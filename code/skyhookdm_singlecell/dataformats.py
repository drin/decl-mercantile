"""
Sub-module that contains code related to serialization and conversion into various data formats.
For example, serialize a GeneExpression object into Arrow format (using RecordBatches).
"""

# core libraries
import sys
import logging

# dependencies
import numpy
import pyarrow

# classes from this package
from skyhookdm_singlecell import skyhook
from skyhookdm_singlecell.datatypes import GeneExpression

# variables from this package
from skyhookdm_singlecell import (__skyhook_version__            ,
                                  __skyhook_data_schema_version__,
                                  __skyhook_data_struct_version__,)

# ------------------------------
# Module-level Variables
debug = True


# ------------------------------
# Classes
class SkyhookSerializerGeneExpression(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)
    # logger.addHandler(logging.StreamHandler(sys.stderr))

    def __init__(self, gene_expression, **kwargs):
        super().__init__(**kwargs)

        self.logger    = self.__class__.logger
        self.gene_expr = gene_expression

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = ''
        self.table_name     = 'gene_expression'

    def arrow_schema(self):
        """
        Create an arrow schema for every column (expression for a cell)

        NOTE: pyarrow data type should match the data type used in skyhook data schema
        """

        return pyarrow.schema((
            #(cell_id, pyarrow.uint16())
            (cell_id, pyarrow.float16())
            for cell_id in self.gene_expr.cells
        ))

    def skyhook_column_schemas(self):
        return [
            skyhook.ColumnSchema(
                cell_ndx                       ,
                skyhook.DataTypes.SDT_UINT16   ,
                skyhook.KeyColumn.NOT_KEY      ,
                skyhook.NullableColumn.NULLABLE,
                cell_id
            )
            for cell_ndx, cell_id in enumerate(self.gene_expr.cells)
        ]

    def skyhook_metadata(self):
        """
        Return formatted Skyhook metadata that will be stored with the arrow schema information.
        """

        # closure to convert list of column schemas to a single data schema
        def to_skyhook_data_schema(*column_schema):
            return ' '.join(map(str, *column_schema))

        # apply closure to get data schema
        skyhook_data_schema = '\n'.join(map(to_skyhook_data_schema, self.skyhook_column_schemas()))

        # return the Skyhook metadata tuple
        return skyhook.SkyhookMetadata(
            __skyhook_version__               ,
            __skyhook_data_schema_version__   ,
            __skyhook_data_struct_version__   ,
            skyhook.FormatTypes.SFT_ARROW     ,
            skyhook_data_schema               ,
            self.db_schema                    ,
            self.table_name                   ,
            self.gene_expr.expression.shape[0]
        )

    def to_arrow_recordbatch(self, data_schema):
        expr_data            = self.gene_expr.expression.toarray().astype(numpy.float16)
        row_count, col_count = self.gene_expr.expression.shape

        return pyarrow.RecordBatch.from_arrays(
             [
                 pyarrow.array(expr_data[:, barcode_ndx])
                 for barcode_ndx in range(col_count)
             ]
            ,schema=data_schema
        )

    def write_data_as_arrow(self, batch_size, batch_offset, path_to_outfile):
        # get the arrow_schema and set the skyhook metadata
        batch_ndx            = 0
        skyhook_schema_arrow = (
            self.arrow_schema()
                #.with_metadata(self.skyhook_metadata().to_bytes())
        )

        self.logger.info('>>> converting data to Arrow format')
        data_recordbatch = self.to_arrow_recordbatch(skyhook_schema_arrow)
        self.logger.info('<<< data format converted')

        self.logger.debug('{} rows in recordbatch'.format(data_recordbatch.num_rows))
        self.logger.debug('batch size: {}'.format(batch_size))
        self.logger.debug('batch offset: {}'.format(batch_offset))

        batch_writer = pyarrow.RecordBatchFileWriter(
            path_to_outfile,
            schema=skyhook_schema_arrow
        )

        batch_writer.write_batch(data_recordbatch)
        # batch_writer.write_batch(data_recordbatch.slice(offset=0, length=batch_size))

        '''
        while batch_ndx + batch_offset < data_recordbatch.num_rows:

            if debug:
                sys.stdout.write('\rwriting batch: {:0d}'.format(batch_ndx))
                sys.stdout.flush()

            batch_writer.write_batch(data_recordbatch.slice(offset=batch_ndx, length=batch_size))
            batch_ndx += batch_offset

        # have to serialize the remainder batch (num_rows % batch_offset)
        if data_recordbatch.num_rows > batch_ndx:
            sys.stdout.write('\rwriting batch: {:0d}\n'.format(batch_ndx))
            sys.stdout.flush()

            batch_writer.write_batch(data_recordbatch.slice(offset=batch_ndx, length=batch_size))
        '''

        self.logger.info('--- batches written')
