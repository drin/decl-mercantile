"""
Sub-module that contains code related to serialization and conversion into various data formats.
For example, serialize a GeneExpression object into Arrow format (using RecordBatches).
"""

# core libraries
import os
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


# Expression Value Representation
# NOTE: these need to match logically, even though they have different representations
numpy_type   = numpy.uint16
arrow_type   = pyarrow.uint16()
skyhook_type = skyhook.DataTypes.SDT_UINT16


# ------------------------------
# Static data conversion functions
def arrow_binary_from_table(arrow_table):
    arrow_buffer  = pyarrow.BufferOutputStream()
    stream_writer = pyarrow.RecordBatchStreamWriter(arrow_buffer, arrow_table.schema)

    for record_batch in arrow_table.to_batches():
        stream_writer.write_batch(record_batch)

    stream_writer.close()

    # TODO: see if to_pybytes() is necessary
    # return arrow_buffer.getvalue().to_pybytes()

    return arrow_buffer.getvalue()


def arrow_table_from_binary(cls, data_blob):
    stream_reader = pyarrow.ipc.open_stream(data_blob)

    deserialized_batches = [
        deserialized_batch
        for deserialized_batch in stream_reader
    ]

    stream_reader.close()

    return pyarrow.Table.from_batches(deserialized_batches, schema=stream_reader.schema)


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
        self.gene_expr.expression.astype(dtype=numpy_type)

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = ''
        self.table_name     = 'gene_expression'

    def arrow_schema_for_cells(self, start_ndx=None, end_ndx=None):
        # NOTE: pyarrow data type should match the data type used in skyhook data schema
        return pyarrow.schema([
            (cell_id, arrow_type)
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
                col_id                         ,
                skyhook_type                   ,
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
                        .astype(dtype=numpy_type),
                    cell_stop_ndx - cell_start_ndx
                )
            ))

        else:
            cell_expr = (
                self.gene_expr
                    .expression
                    .transpose()
                    .toarray()
                    .astype(dtype=numpy_type)
            )

        return pyarrow.RecordBatch.from_arrays(cell_expr, schema=data_schema)

    def to_arrow_table(self, data_schema):
        # convert matrix or ndarray into a list of 1-d arrays
        cell_expression_arrays = list(
            # apply `ravel()` to each list from `hsplit()` to yield list of <array_like>
            map(
                numpy.ravel,
                # `hsplit()` returns a list of ndarrays
                numpy.hsplit(
                    self.gene_expr.expression.astype(dtype=numpy_type),
                    self.gene_expr.expression.shape[1]
                )
            )
        )

        return pyarrow.Table.from_arrays(cell_expression_arrays, schema=data_schema)

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
                    .astype(dtype=numpy_type)
            )

        elif isinstance(self.gene_expr.expression, numpy.matrix):
            cell_expr = (
                self.gene_expr
                    .expression[:, start:end]
                    .getA()
                    .astype(dtype=numpy_type)
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


class SkyhookFlatbufferMeta(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_binary_flatbuffer(cls, flatbuffer_binary, offset=0):
        return cls(FB_Meta.FB_Meta.GetRootAsFB_Meta(flatbuffer_binary, offset))

    @classmethod
    def binary_from_arrow_table(cls, binary_arrow_table):
        # initialize a flatbuffer builder with a count of expected contiguous bytes needed,
        # accommodating each fixed-size field of the FB_Meta flatbuffer
        partial_byte_count = (4 + 8 + 4 + 8 + 8 + 4)
        builder = flatbuffers.Builder(partial_byte_count + len(binary_arrow_table))

        # add the serialized data first (build flatbuffer from back to front)
        #wrapped_data_blob = builder.CreateByteVector(binary_arrow_table.to_pybytes())
        wrapped_data_blob = builder.CreateByteVector(binary_arrow_table)

        # construct the remaining flatbuffer structure
        FB_Meta.FB_MetaStart(builder)

        FB_Meta.FB_MetaAddBlobFormat(builder, skyhook.FormatTypes.SFT_ARROW)
        FB_Meta.FB_MetaAddBlobData(builder, wrapped_data_blob)
        FB_Meta.FB_MetaAddBlobSize(builder, len(binary_arrow_table))
        FB_Meta.FB_MetaAddBlobDeleted(builder, False)
        FB_Meta.FB_MetaAddBlobOrigOff(builder, 0)
        FB_Meta.FB_MetaAddBlobOrigLen(builder, len(binary_arrow_table))
        FB_Meta.FB_MetaAddBlobCompression(builder, 0)

        builder.Finish(FB_Meta.FB_MetaEnd(builder))

        # return the finished binary blob representing FBMeta(<Arrow binary>)
        return builder.Output()

    def __init__(self, flatbuffer_obj, **kwargs):
        super().__init__(**kwargs)

        self.fb_obj = flatbuffer_obj

    def get_data_format(self):
        return self.fb_obj.BlobFormat() or None

    def get_size(self):
        return self.fb_obj.BlobDataLength() or None

    def get_data_as_arrow(self):
        data_blob = self.fb_obj.BlobDataAsNumpy()

        if type(data_blob) is int and data_blob == 0: return None

        stream_reader      = pyarrow.ipc.open_stream(data_blob.tobytes())
        deserialized_table = pyarrow.Table.from_batches(
            list(stream_reader),
            schema=stream_reader.schema
        )

        return deserialized_table


class SkyhookFileReader(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def read_gene_expr_partitions_to_arrow_table(cls, path_to_directory, file_ext='.arrow'):
        cls.logger.info('>>> reading data partitions from directory of arrow files')

        # initialize these variables, just in case there are no files in path_to_directory
        table_partitions = []

        for partition_file in os.listdir(path_to_directory):
            if not partition_file.endswith(file_ext): continue

            with open(os.path.join(path_to_directory, partition_file), 'rb') as input_handle:
                table_partitions.append(arrow_table_from_binary(input_handle.read()))

        cls.logger.info(f'<<< read {len(table_partitions)} partitions')
        return pyarrow.concat_tables(table_partitions)

    @classmethod
    def read_gene_expr_to_arrow_table(cls, path_to_infile):
        cls.logger.info('>>> reading gene expression from single arrow file')

        with open(path_to_infile, 'rb') as input_handle:
            gene_expr_table = arrow_table_from_binary(input_handle.read())

        cls.logger.info('<<< data read into arrow table')

        return gene_expr_table

    @classmethod
    def read_gene_expr_to_arrow_binary(cls, path_to_infile):
        cls.logger.info('>>> reading gene expression from single arrow file')

        with open(path_to_infile, 'rb') as input_handle:
            gene_expr_table = input_handle.read()

        cls.logger.info('<<< data read as binary')

        return gene_expr_table


class SkyhookFileWriter(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def write_gene_expr_partitions_flatbuffer(cls, gene_expr, output_dir, batch_size=100):
        # tbl_part is table partition
        for ndx, tbl_part in gene_expr.cell_expression_batched(batch_size):
            # Generate path to binary file based on table partition metadata
            path_to_blob = os.path.join(
                output_dir,
                '{}-{}-{}.skyhook'.format(ndx, tbl_part.num_columns, tbl_part.column_names[0])
            )

            # Serialize flatbuffer structure to disk
            fb_meta_wrapped_binary = SkyhookFlatbufferMeta.binary_from_arrow_table(
                arrow_binary_from_table(tbl_part)
            )
            with open(path_to_blob, 'wb') as blob_handle:
                blob_handle.write(fb_meta_wrapped_binary)

    @classmethod
    def write_gene_expr_partitions_arrow(cls, gene_expr, output_dir, batch_size=100):
        cls.logger.info('>>> writing data partitions into directory of arrow files')

        # tbl_part is table partition
        for ndx, tbl_part in gene_expr.cell_expression_batched(batch_size=batch_size):
            # Generate path to binary file based on table partition metadata
            path_to_blob = os.path.join(
                output_dir,
                '{}-{}-{}.arrow'.format(ndx, tbl_part.num_columns, tbl_part.column_names[0])
            )

            # Serialize flatbuffer structure to disk
            arrow_binary_data = arrow_binary_from_table(tbl_part)
            with open(path_to_blob, 'wb') as blob_handle:
                blob_handle.write(arrow_binary_data)

        cls.logger.info('<<< data written')

    @classmethod
    def write_gene_expr_arrow(cls, gene_expr, path_to_outfile):
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

    @classmethod
    def write_gene_expr_partitions_parquet(cls, gene_expr, path_to_directory, batch_size=1000):
        cls.logger.info('>>> writing data partitions into directory of parquet files')

        # tbl_part is table partition
        for ndx, tbl_part in gene_expr.cell_expression_batched(batch_size=batch_size):
            path_to_parquet_file = os.path.join(path_to_directory, '{}-{}-{}.parquet'.format(
                ndx, tbl_part.num_columns, tbl_part.column_names[0]
            ))

            pyarrow.parquet.write_table(tbl_part, path_to_parquet_file)

        cls.logger.info('<<< data written')

    @classmethod
    def write_gene_expr_parquet(cls, gene_expr, path_to_outfile):
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
