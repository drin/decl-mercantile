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

from pyarrow import parquet

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

    def __init__(self, gene_expression, **kwargs):
        super().__init__(**kwargs)

        self.logger    = self.__class__.logger
        self.gene_expr = gene_expression

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = ''
        self.table_name     = 'gene_expression'

    def schema_for_cell_expression(self, col_id, cell_id):
        return skyhook.ColumnSchema(
            col_id,
            skyhook.DataTypes.SDT_UINT16   ,
            skyhook.KeyColumn.NOT_KEY      ,
            skyhook.NullableColumn.NULLABLE,
            cell_id
        )

    def schema_for_gene_expression(self, col_id, gene_name):
        return skyhook.ColumnSchema(
            col_id,
            skyhook.DataTypes.SDT_UINT16   ,
            skyhook.KeyColumn.NOT_KEY      ,
            skyhook.NullableColumn.NULLABLE,
            gene_name
        )

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

        # apply closure to get data schema
        skyhook_data_schema = '\n'.join(map(str, self.skyhook_column_schemas()))

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

    def skyhook_metadata_for_cell_expression(self, col_id, cell_id):
        """
        Return formatted Skyhook metadata that will be stored with the arrow schema information.
        """

        cell_expr_schema = self.schema_for_cell_expression(col_id, cell_id)

        # return the Skyhook metadata tuple
        return skyhook.SkyhookMetadata(
            __skyhook_version__               ,
            __skyhook_data_schema_version__   ,
            __skyhook_data_struct_version__   ,
            skyhook.FormatTypes.SFT_ARROW     ,
            str(cell_expr_schema)             ,
            self.db_schema                    ,
            self.table_name                   ,
            self.gene_expr.expression.shape[0]
        )

    def arrow_schema(self):
        """
        Create an arrow schema for every column (expression for a cell)

        NOTE: pyarrow data type should match the data type used in skyhook data schema
        """

        return pyarrow.schema((
            (cell_id, pyarrow.int8())
            for cell_id in self.gene_expr.cells
        ))

    def arrow_schema_for_cell(self, cell_id):
        return pyarrow.schema([ (cell_id, pyarrow.int8()) ])

    def to_arrow_recordbatch(self, data_schema):
        """
        Arrow format assumes that the columns are in the first dimension, not the second dimension.

        So we pass the expression numpy array to pyarrow as a transposed matrix (cells are rows,
        genes are columns). The schema is a column-based schema (cells).
        """

        return pyarrow.RecordBatch.from_arrays(self.gene_expr.expression.T, schema=data_schema)

    def to_arrow_table(self, data_schema):
        return pyarrow.Table.from_arrays(self.gene_expr.expression.T, schema=data_schema)

    def serialize_skyhook_metadata(self):
        return self.skyhook_metadata().to_byte_coercible()

    def deserialize_skyhook_metadata(self, serialized_metadata):
        return skyhook.SkyhookMetadata(**serialized_metadata)

    def arrow_serialization_context(self):
        context = pyarrow.SerializationContext()

        context.register_type(
            skyhook.SkyhookMetadata,
            'SkyhookMetadata',
            custom_serializer=self.serialize_skyhook_metadata,
            custom_deserializer=self.deserialize_skyhook_metadata
        )

        return context

    def write_data_as_arrow(self, batch_size, batch_offset, path_to_outfile):
        # get the arrow_schema and set the skyhook metadata
        batch_ndx            = 0
        skyhook_schema_arrow = (
            self.arrow_schema()
                .with_metadata(self.serialize_skyhook_metadata())
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

        self.logger.info('--- batches written')

    def write_data_as_parquet(self, path_to_outfile):
        data_schema = (
            self.arrow_schema()
                .with_metadata(self.serialize_skyhook_metadata())
        )

        pyarrow.parquet.write_table(self.to_arrow_table(data_schema), path_to_outfile)

    def write_cells_as_parquet(self, path_to_directory):
        cell_expr = self.gene_expr.expression.T

        for cell_ndx, cell_id in enumerate(self.gene_expr.cells):
            cell_expr_schema = (
                self.arrow_schema_for_cell(cell_id)
                    .with_metadata(
                        self.skyhook_metadata_for_cell_expression(cell_ndx, cell_id)
                            .to_byte_coercible()
                     )
            )

            pyarrow.parquet.write_table(
                pyarrow.Table.from_arrays(
                    (cell_expr[cell_ndx],), schema=cell_expr_schema
                ),
                os.path.join(path_to_directory, f'{cell_expr}.parquet')
            )
