"""
Sub-module that contains code related to serialization and conversion into various data formats.
For example, serialize a GeneExpression object into Arrow format (using RecordBatches).
"""

# core libraries
import os
import sys
import logging

# dependencies
import numpy
import scipy
import pyarrow
import pyarrow.parquet

# classes from this package
from skyhookdm_singlecell import skyhook

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

        # hardcoded for now
        self.gene_expr.expression.astype(dtype=numpy.int8)

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = ''
        self.table_name     = 'gene_expression'

    def schema_for_cell_expression(self, col_ids, cell_ids):
        return [
            skyhook.ColumnSchema(
                col_id,
                skyhook.DataTypes.SDT_UINT16   ,
                skyhook.KeyColumn.NOT_KEY      ,
                skyhook.NullableColumn.NULLABLE,
                cell_id
            )

            for col_id, cell_id in zip(col_ids, cell_ids)
        ]

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

    def metadata_by_cells(self, col_ids, cell_ids):
        """
        Return formatted Skyhook metadata that will be stored with the arrow schema information.
        """

        cell_expr_schema = '\n'.join(map(str, self.schema_for_cell_expression(col_ids, cell_ids)))

        # return the Skyhook metadata tuple
        return skyhook.SkyhookMetadata(
            __skyhook_version__               ,
            __skyhook_data_schema_version__   ,
            __skyhook_data_struct_version__   ,
            skyhook.FormatTypes.SFT_ARROW     ,
            cell_expr_schema                  ,
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

    def schema_for_cells(self, cell_ids):
        return pyarrow.schema([
            (cell_id, pyarrow.int8())
            for cell_id in cell_ids
        ])

    def cells_to_arrow_recordbatch(self, data_schema, cell_start_ndx=None, cell_stop_ndx=None):
        """
        Arrow format assumes that the columns are in the first dimension, not the second dimension.

        So we pass the expression numpy array to pyarrow as a transposed matrix (cells are rows,
        genes are columns). The schema is a column-based schema (cells).
        """

        if cell_start_ndx is not None and cell_stop_ndx is not None:
            cell_expr = list(map(
                numpy.ravel,
                numpy.hsplit(
                    self.gene_expr
                        .expression[:, cell_start_ndx:cell_stop_ndx]
                        .getA()
                        .astype(dtype=numpy.int8),
                    cell_stop_ndx - cell_start_ndx
                )
            ))

        else:
            cell_expr = (
                self.gene_expr
                    .expression
                    .transpose()
                    .toarray()
                    .astype(dtype=numpy.int8)
            )

        return pyarrow.RecordBatch.from_arrays(cell_expr, schema=data_schema)

    def to_arrow_table(self, data_schema):
        cell_expression_arrays = [
            self.gene_expr.expression.toarray()[:, cell_ndx]
            for cell_ndx in range(self.gene_expr.expression.shape[1])
        ]

        return pyarrow.Table.from_arrays(cell_expression_arrays, schema=data_schema)

    '''
    Deprecated for now
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
    '''

    def serialize_as_sft_arrow(self, arrow_table):
        arrow_buffer  = pyarrow.BufferOutputStream()
        stream_writer = pyarrow.RecordBatchStreamWriter(arrow_buffer, arrow_table.schema)

        for record_batch in arrow_table.to_batches():
            stream_writer.write_batch(record_batch)

        stream_writer.close()

        return arrow_buffer.getvalue()

    def deserialize_from_sft_arrow(self, data_blob):
        stream_reader = pyarrow.ipc.open_stream(data_blob)

        deserialized_batches = [
            deserialized_batch
            for deserialized_batch in stream_reader
        ]

        stream_reader.close()

        return deserialized_batches

    def cell_expression_table(self, start=0, end=None):
        if end is None:
            end = self.gene_expr.cells.shape[0]

        cell_indices  = range(start, end)
        cell_ids      = self.gene_expr.cells[start:end]

        # ------------------------------
        # Construct schema
        self.logger.info(f'>>> constructing schema for partition ({start}:{end})')

        skyhook_metadata = self.metadata_by_cells(cell_indices, cell_ids).to_byte_coercible()
        cell_expr_schema = (
            self.schema_for_cells(cell_ids)
                .with_metadata(skyhook_metadata)
        )

        self.logger.info('<<< partition schema constructed')

        # ------------------------------
        # slice of cell expression (columns) as arrays
        self.logger.info('>>> extracting cell expression slice')

        cell_expr = None

        if isinstance(self.gene_expr.expression, scipy.sparse.coo_matrix):
            cell_expr = (
                self.gene_expr
                    .expression
                    .toarray()[:, start:end]
                    .astype(dtype=numpy.int8)
            )

        elif isinstance(self.gene_expr.expression, numpy.matrix):
            cell_expr = (
                self.gene_expr
                    .expression[:, start:end]
                    .getA()
                    .astype(dtype=numpy.int8)
            )

        self.logger.info('<<< cell expression slice extracted')

        return pyarrow.Table.from_arrays(
            list(map(numpy.ravel, numpy.hsplit(cell_expr, cell_ids.shape[0]))),
            schema=cell_expr_schema
        )

    def batched_indices(self, element_count, batch_size):
        for batch_id, batch_start in enumerate(range(0, element_count, batch_size), start=1):

            batch_end = batch_id * batch_size
            if batch_end > element_count:
                batch_end  = element_count

            yield batch_id, batch_start, batch_end

    def cell_expression_batched(self, batch_size=1000, cell_count=None):
        if cell_count is None:
            cell_count = self.gene_expr.cells.shape[0]

        for batch_id, start, end in self.batched_indices(cell_count, batch_size):
            yield batch_id, self.cell_expression_table(start=start, end=end)

    def write_cells_as_arrow(self, path_to_directory, batch_size=1000):
        for partition_ndx, cell_expr_table in self.cell_expression_batched(batch_size=batch_size):
            path_to_arrow_file = os.path.join(path_to_directory, '{}-{}-{}.arrow'.format(
                partition_ndx, cell_expr_table.num_columns, cell_expr_table.column_names[0]
            ))

            with pyarrow.OSFile(path_to_arrow_file, 'wb') as arrow_handle:
                batch_writer = pyarrow.RecordBatchFileWriter(arrow_handle, cell_expr_table.schema)

                for record_batch in cell_expr_table.to_batches():
                    batch_writer.write_batch(record_batch)

    def write_cells_as_parquet(self, path_to_directory, batch_size=1000):
        for partition_ndx, cell_expr_table in self.cell_expression_batched(batch_size=batch_size):
            path_to_parquet_file = os.path.join(path_to_directory, '{}-{}-{}.parquet'.format(
                partition_ndx, cell_expr_table.num_columns, cell_expr_table.column_names[0]
            ))

            pyarrow.parquet.write_table(cell_expr_table, path_to_parquet_file)

    def write_data_as_arrow(self, path_to_outfile):
        data_schema = (
            self.arrow_schema()
                .with_metadata(self.serialize_skyhook_metadata())
        )

        self.logger.info('>>> writing data in single arrow file')

        with open(path_to_outfile, 'wb') as arrow_handle:
                batch_writer = pyarrow.RecordBatchFileWriter(arrow_handle, data_schema)

                for record_batch in self.to_arrow_table(data_schema).to_batches():
                    batch_writer.write_batch(record_batch)

        self.logger.info('<<< data written')

    def write_data_as_parquet(self, path_to_outfile):
        data_schema = (
            self.arrow_schema()
                .with_metadata(self.serialize_skyhook_metadata())
        )

        self.logger.info('>>> writing data in single parquet file')

        pyarrow.parquet.write_table(self.to_arrow_table(data_schema), path_to_outfile)

        self.logger.info('<<< data written')
