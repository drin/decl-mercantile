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
import flatbuffers

# modules from this package
from skyhookdm_singlecell import skyhook
from skyhookdm_singlecell.Tables import Table, Record, FB_Meta

# variables from this package
from skyhookdm_singlecell import (__skyhook_version__            ,
                                  __skyhook_data_schema_version__,
                                  __skyhook_data_struct_version__,)

# ------------------------------
# Module-level Variables
debug = True


# ------------------------------
# Classes
class SkyhookGeneExpression(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    # TODO: make the data type to serialize and format to available module wide and/or to
    # everything that needs to do the conversion. Hopefully parameterize it too

    def __init__(self, gene_expression, **kwargs):
        super().__init__(**kwargs)

        self.logger    = self.__class__.logger
        self.gene_expr = gene_expression

        # hardcoded for now
        self.gene_expr.expression.astype(dtype=numpy.uint8)

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = ''
        self.table_name     = 'gene_expression'

    def arrow_schema_for_cells(self, start_ndx=None, end_ndx=None):
        # NOTE: pyarrow data type should match the data type used in skyhook data schema
        return pyarrow.schema([
            (cell_id, pyarrow.uint8())
            for cell_id in self.gene_expr.cells[start_ndx:end_ndx]
        ])

    def schema_by_cell(self, start_ndx=0, end_ndx=None):
        """
        Create the skyhook schema description for a subset of the cells in this gene expression
        matrix. Note, the `DataType` field should be kept in-sync with the datatype used in
        serialization (see `arrow_schema_for_cells`).
        """

        # NOTE: this iterator can only be walked once
        cell_iterator = enumerate(
            self.gene_expr.cells[start_ndx:end_ndx],
            start=start_ndx
        )

        return [
            skyhook.ColumnSchema(
                col_id,
                skyhook.DataTypes.SDT_UINT8    ,
                skyhook.KeyColumn.NOT_KEY      ,
                skyhook.NullableColumn.NULLABLE,
                cell_id
            )

            for col_id, cell_id in cell_iterator
        ]

    def skyhook_metadata_by_cell(self, start_ndx=0, end_ndx=None):
        """
        Return formatted Skyhook metadata that will be stored with the arrow schema information.
        """

        cell_expr_schema = '\n'.join(map(
            str,
            self.schema_by_cell(start_ndx=start_ndx, end_ndx=end_ndx)
        ))

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
                        .astype(dtype=numpy.uint8),
                    cell_stop_ndx - cell_start_ndx
                )
            ))

        else:
            cell_expr = (
                self.gene_expr
                    .expression
                    .transpose()
                    .toarray()
                    .astype(dtype=numpy.uint8)
            )

        return pyarrow.RecordBatch.from_arrays(cell_expr, schema=data_schema)

    def to_arrow_table(self, data_schema):
        cell_expression_arrays = [
            self.gene_expr.expression.toarray()[:, cell_ndx]
            for cell_ndx in range(self.gene_expr.expression.shape[1])
        ]

        return pyarrow.Table.from_arrays(cell_expression_arrays, schema=data_schema)

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

        # ------------------------------
        # Construct schema
        self.logger.info(f'>>> constructing schema for partition ({start}:{end})')

        skyhook_metadata = self.skyhook_metadata_by_cell(start, end).to_byte_coercible()
        cell_expr_schema = (
            self.arrow_schema_for_cells(start_ndx=start, end_ndx=end)
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
                    .astype(dtype=numpy.uint8)
            )

        elif isinstance(self.gene_expr.expression, numpy.matrix):
            cell_expr = (
                self.gene_expr
                    .expression[:, start:end]
                    .getA()
                    .astype(dtype=numpy.uint8)
            )

        self.logger.info('<<< cell expression slice extracted')

        return pyarrow.Table.from_arrays(
            list(map(numpy.ravel, numpy.hsplit(cell_expr, cell_expr.shape[1]))),
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


class GeneExpressionArrowWriter(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def write_cells_as_arrow(cls, gene_expr, path_to_directory, batch_size=1000):
        cls.logger.info('>>> writing data partitions into directory of arrow files')

        for ndx, table_partition in gene_expr.cell_expression_batched(batch_size=batch_size):
            arrow_buffer = pyarrow.BufferOutputStream()

            # write record batches to arrow_buffer
            batch_writer = pyarrow.RecordBatchStreamWriter(arrow_buffer, table_partition.schema)
            for record_batch in table_partition.to_batches():
                batch_writer.write_batch(record_batch)
            batch_writer.close()

            # write contents of arrow_buffer to binary file
            path_to_arrow_file = os.path.join(path_to_directory, '{}-{}-{}.arrow'.format(
                ndx, table_partition.num_columns, table_partition.column_names[0]
            ))
            with open(path_to_arrow_file, 'wb') as output_handle:
                output_handle.write(arrow_buffer.getvalue().to_pybytes())

        cls.logger.info('<<< data written')

    @classmethod
    def write_data_as_arrow(cls, gene_expr, path_to_outfile):
        data_schema = (
            gene_expr.arrow_schema_for_cells()
                     .with_metadata(
                          gene_expr.skyhook_metadata_by_cell()
                                   .to_byte_coercible()
                      )
        )

        cls.logger.info('>>> writing data in single arrow file')

        with open(path_to_outfile, 'wb') as arrow_handle:
                batch_writer = pyarrow.RecordBatchFileWriter(arrow_handle, data_schema)

                for record_batch in gene_expr.to_arrow_table(data_schema).to_batches():
                    batch_writer.write_batch(record_batch)

        cls.logger.info('<<< data written')


class GeneExpressionParquetWriter(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def write_data_as_parquet(cls, gene_expr, path_to_outfile):
        data_schema = (
            gene_expr.arrow_schema_for_cells()
                     .with_metadata(
                          gene_expr.skyhook_metadata_by_cell()
                                   .to_byte_coercible()
                      )
        )

        cls.logger.info('>>> writing data in single parquet file')

        pyarrow.parquet.write_table(gene_expr.to_arrow_table(data_schema), path_to_outfile)

        cls.logger.info('<<< data written')

    @classmethod
    def write_cells_as_parquet(cls, gene_expr, path_to_directory, batch_size=1000):
        cls.logger.info('>>> writing data partitions into directory of parquet files')

        for ndx, table_partition in gene_expr.cell_expression_batched(batch_size=batch_size):
            path_to_parquet_file = os.path.join(path_to_directory, '{}-{}-{}.parquet'.format(
                ndx, table_partition.num_columns, table_partition.column_names[0]
            ))

            pyarrow.parquet.write_table(table_partition, path_to_parquet_file)

        cls.logger.info('<<< data written')


class SkyhookFlatbufferTable(object):
    """
    Flatbuffer methods return 0 on failure.
    """

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_binary(cls, flatbuffer_binary, offset=0):
        return cls(Table.Table.GetRootAsTable(flatbuffer_binary, offset))

    def __init__(self, flatbuffer_obj, **kwargs):
        super().__init__(**kwargs)

        self.fb_obj = flatbuffer_obj

    def get_data_format(self):
        return skyhook.FormatTypes.__enum_names__[self.fb_obj.DataFormatType()]

    def is_arrow_type(self):
        return self.fb_obj.DataFormatType() == skyhook.FormatTypes.SFT_ARROW

    def get_schema(self):
        return self.fb_obj.DataSchema()

    def get_table_name(self):
        return self.fb_obj.TableName()

    def get_data_as_arrow(self):
        pass

    def get_row_count(self):
        return self.fb_obj.Nrows()


class SkyhookFlatbufferMeta(object):
    """
    Flatbuffer methods return 0 on failure.
    """

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_binary(cls, flatbuffer_binary, offset=0):
        return cls(FB_Meta.FB_Meta.GetRootAsFB_Meta(flatbuffer_binary, offset))

    def __init__(self, flatbuffer_obj, **kwargs):
        super().__init__(**kwargs)

        self.fb_obj = flatbuffer_obj

    def get_data_format(self):
        return self.fb_obj.BlobFormat() or None

    def get_size(self):
        return self.fb_obj.BlobDataLength() or None

    def get_data_as_arrow(self):
        data_blob = self.fb_obj.BlobDataAsNumpy()

        if not data_blob: return None

        stream_reader      = pyarrow.ipc.open_stream(data_blob.tobytes())
        deserialized_table = pyarrow.Table.from_batches(
            list(stream_reader),
            schema=stream_reader.schema
        )

        return deserialized_table

    '''
    def get_data_as_flex(self):
        data_blob = self.fb_obj.BlobDataAsNumpy()

        if not data_blob: return None

        return SkyhookFlatbufferTable.from_binary(data_blob.tobytes())
    '''


class SkyhookFlatbufferWriter(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def tabular_arrow_buffer(self, gene_expr, output_dir, batch_size=500):
        """
        Write object, in-order
        """

        # byte counts of fixed-size, or known apriori, fields
        partial_byte_count = (
            4 + 4 + 4 + 4 +
            len(gene_expr.table_name.encode('utf-8')) +
            8 + 4 + 8 + 8 + 4
        )

        # for each partition of the gene expression (arrow table)
        for ndx, table_partition in gene_expr.cell_expression_batched(batch_size):
            table_binary = gene_expr.serialize_as_sft_arrow(table_partition)

            # byte counts of fields that needed to be computed
            table_byte_count = (
                len(table_partition.schema.metadata['METADATA_DATA_SCHEMA'.encode('utf-8')]) +
                (table_partition.num_columns << 4) +
                table_binary.size
            )

            builder = flatbuffers.Builder(partial_byte_count + table_byte_count)

            # we want to serialize the data first (build vector from back to front)
            wrapped_data_blob = builder.CreateByteVector(table_binary.to_pybytes())

            FB_Meta.FB_MetaStart(builder)

            FB_Meta.FB_MetaAddBlobFormat(builder, 1)
            FB_Meta.FB_MetaAddBlobData(builder, wrapped_data_blob)
            FB_Meta.FB_MetaAddBlobSize(builder, table_binary.size)
            FB_Meta.FB_MetaAddBlobDeleted(builder, False)
            FB_Meta.FB_MetaAddBlobOrigOff(builder, 0)
            FB_Meta.FB_MetaAddBlobOrigLen(builder, table_binary.size)
            FB_Meta.FB_MetaAddBlobCompression(builder, 0)

            builder.Finish(FB_Meta.FB_MetaEnd(builder))

            return builder.Output()

    @classmethod
    def table_record_arrow_buffer(self, row_id, column_count, arrow_buffer, initial_byte_count=0):

        partial_byte_count = 8

        if not initial_byte_count:
            # ID + 8 bytes per field + bytes of data
            initial_byte_count = 8 + (column_count << 3)  + arrow_buffer.size

        # for each partition of the gene expression (arrow table)
        for ndx, table_partition in gene_expr.cell_expression_batched(batch_size):
            table_binary = gene_expr.serialize_as_sft_arrow(table_partition)

        builder = flatbuffers.Builder(initial_byte_count)

        Record.RecordStart(builder)

        # RID is a uint64, and so it's 8 bytes
        Record.RecordAddRID(row_id.to_bytes(8, byteorder='big'))

        # Record.RecordAddNullbits

        builder.Finish(Record.RecordEnd(builder))

        return builder.Output()

    @classmethod
    def wrapped_arrow_buffer(self, gene_expr, output_dir, batch_size=500):
        """
        Write object, in-order
        """

        # accommodate each fixed-size field of the FB_Meta flatbuffer
        partial_byte_count = (4 + 8 + 4 + 8 + 8 + 4)

        # for each partition of the gene expression (arrow table)
        for ndx, table_partition in gene_expr.cell_expression_batched(batch_size):
            table_binary = gene_expr.serialize_as_sft_arrow(table_partition)

            # construct the flatbuffer structure to contain the binary data
            builder = flatbuffers.Builder(partial_byte_count + table_binary.size)

            # add the serialized data first (build flatbuffer from back to front)
            wrapped_data_blob = builder.CreateByteVector(table_binary.to_pybytes())

            FB_MetaStart(builder)

            FB_MetaAddBlobFormat(builder, 1)
            FB_MetaAddBlobData(builder, wrapped_data_blob)
            FB_MetaAddBlobSize(builder, table_binary.size)
            FB_MetaAddBlobDeleted(builder, False)
            FB_MetaAddBlobOrigOff(builder, 0)
            FB_MetaAddBlobOrigLen(builder, table_binary.size)
            FB_MetaAddBlobCompression(builder, 0)

            builder.Finish(FB_MetaEnd(builder))

            # Serialize flatbuffer structure to disk
            path_to_blob = os.path.join(
                output_dir,
                '{}-{}-{}.skyhook-blob'.format(
                    ndx, table_partition.num_columns, table_partition.column_names[0]
                )
            )

            with open(path_to_blob, 'wb') as blob_handle:
                blob_handle.write(builder.Output())
