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
from skyhookdm_singlecell.Tables import FB_Meta

# functions
from skyhookdm_singlecell.util import batched_indices

# variables from this package
from skyhookdm_singlecell import (__skyhook_version__            ,
                                  __skyhook_data_schema_version__,
                                  __skyhook_data_struct_version__,)

# ------------------------------
# Module-level Variables
debug = True


# ------------------------------
# Static data conversion functions
def arrow_binary_from_table(arrow_table):
    arrow_buffer  = pyarrow.BufferOutputStream()
    stream_writer = pyarrow.RecordBatchStreamWriter(arrow_buffer, arrow_table.schema)

    for record_batch in arrow_table.to_batches():
        stream_writer.write_batch(record_batch)

    stream_writer.close()

    return arrow_buffer.getvalue()


def arrow_batches_from_binary(data_blob):
    stream_reader = pyarrow.ipc.open_stream(data_blob)

    return (
        stream_reader.schema,
        [deserialized_batch for deserialized_batch in stream_reader]
    )


def arrow_table_from_binary(data_blob):
    blob_schema, blob_batches = arrow_batches_from_binary(data_blob)

    return pyarrow.Table.from_batches(blob_batches, schema=blob_schema)


# ------------------------------
# Classes
class SkyhookDataWrapper(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    def __init__(self, table_name, domain_dataset, dataset_columns, db_schema='public',
                 type_for_numpy=numpy.uint16, type_for_arrow=pyarrow.uint16(),
                 type_for_skyhook=skyhook.DataTypes.SDT_UINT16, **kwargs):

        super().__init__(**kwargs)

        # maintain type information for data conversions
        self.numpy_type     = type_for_numpy
        self.arrow_type     = type_for_arrow
        self.skyhook_type   = type_for_skyhook

        self.logger         = self.__class__.logger
        self.domain_data    = domain_dataset.astype(dtype=self.numpy_type, copy=False)
        self.dataset_cols   = dataset_columns

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = db_schema
        self.table_name     = table_name

    def table_schema(self, from_col_ndx=None, to_col_ndx=None):
        return pyarrow.schema([
            (column_name, self.arrow_type)
            for column_name in self.dataset_cols[from_col_ndx:to_col_ndx]
        ])

    def skyhook_data_schema(self, from_col_ndx=0, to_col_ndx=None):
        col_iterator = enumerate(
            self.dataset_cols[from_col_ndx:to_col_ndx],
            start=from_col_ndx
        )

        return [
            skyhook.ColumnSchema(
                column_ndx                     ,
                self.skyhook_type              ,
                skyhook.KeyColumn.NOT_KEY      ,
                skyhook.NullableColumn.NULLABLE,
                column_name
            )

            for column_ndx, column_name in col_iterator
        ]

    def schema_skyhook_metadata(self, from_col_ndx=0, to_col_ndx=None, col_delim=';'):
        """
        Return formatted Skyhook metadata that will be stored with the arrow schema information.
        """

        # data schema is '\n' or ';' delimited information for each column
        data_schema_as_str = col_delim.join(map(
            str,
            self.skyhook_data_schema(from_col_ndx, to_col_ndx)
        ))

        # return the Skyhook metadata tuple
        return skyhook.SkyhookMetadata(
            __skyhook_version__               ,
            __skyhook_data_schema_version__   ,
            __skyhook_data_struct_version__   ,
            skyhook.FormatTypes.SFT_ARROW     ,
            data_schema_as_str                ,
            self.db_schema                    ,
            self.table_name                   ,
            self.domain_data.shape[0]
        )

    def as_arrow_recordbatches(self, data_schema, from_col_ndx=None, to_col_ndx=None):
        """
        Arrow format assumes that the columns are in the first dimension, not the second dimension.

        So we pass the expression numpy array to pyarrow as a transposed matrix (cells are rows,
        genes are columns). The schema is a column-based schema (cells).
        """

        if from_col_ndx is not None and to_col_ndx is not None:
            data_recordbatches = list(map(
                numpy.ravel,
                numpy.hsplit(
                    self.domain_data[:, from_col_ndx:to_col_ndx].getA(),
                    to_col_ndx - from_col_ndx
                )
            ))

        else:
            data_recordbatches = self.domain_data.transpose().toarray()

        return pyarrow.RecordBatch.from_arrays(data_recordbatches, schema=data_schema)

    def as_arrow_table(self, data_schema):
        """
        Creates an Arrow Table with the given data schema.

        Using `hsplit` and `ravel`, the domain_data matrix (or ndarray) is converted into a list of
        1-d arrays. First, `hsplit()` returns a list of ndarrays; second, `ravel()` is applied to
        each list from `hsplit()` to list of <array_like> objects.
        """

        # Converts a matrix (or ndarray) into a list of 1-d arrays.
        column_data_arrays = list(
            map(numpy.ravel, numpy.hsplit(self.domain_data, self.domain_data.shape[1]))
        )

        return pyarrow.Table.from_arrays(column_data_arrays, schema=data_schema)

    def as_partitioned_arrow_table(self, start=0, end=None):
        if end is None: end = self.dataset_cols.shape[0]

        # ------------------------------
        # Construct schema for table partition
        self.logger.info(f'>>> constructing schema for partition ({start}:{end})')

        table_schema_and_meta = (
            self.table_schema(from_col_ndx=start, to_col_ndx=end)
                .with_metadata(
                     self.schema_skyhook_metadata(start, end)
                         .to_byte_coercible()
                 )
        )

        self.logger.info('<<< partition schema constructed')

        # ------------------------------
        # Get slice of domain data as arrays
        self.logger.info('>>> extracting cell expression slice')

        domain_data_as_array = None

        if isinstance(self.domain_data, scipy.sparse.coo_matrix):
            domain_data_as_array = self.domain_data.toarray()[:, start:end]

        elif isinstance(self.domain_data, numpy.matrix):
            domain_data_as_array = self.domain_data[:, start:end].getA()

        elif isinstance(self.domain_data, numpy.ndarray):
            domain_data_as_array = self.domain_data

        self.logger.info('<<< cell expression slice extracted')

        # Create and return arrow table
        return pyarrow.Table.from_arrays(
            list(map(
                numpy.ravel,
                numpy.hsplit(domain_data_as_array, domain_data_as_array.shape[1])
            )),
            schema=table_schema_and_meta
        )

    def batched_table_partitions(self, batch_size=1000, column_count=None):
        if column_count is None: column_count = self.dataset_cols.shape[0]

        for batch_id, start, end in batched_indices(column_count, batch_size):
            yield batch_id, self.as_partitioned_arrow_table(start=start, end=end)


class SkyhookFlatbufferMeta(object):

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_binary_flatbuffer(cls, flatbuffer_binary, offset=0):
        return cls(FB_Meta.FB_Meta.GetRootAsFB_Meta(flatbuffer_binary, offset))

    @classmethod
    def binary_from_arrow_binary(cls, binary_arrow_table):
        # initialize a flatbuffer builder with a count of expected contiguous bytes needed,
        # accommodating each fixed-size field of the FB_Meta flatbuffer
        partial_byte_count = (4 + 8 + 4 + 8 + 8 + 4)
        builder = flatbuffers.Builder(partial_byte_count + len(binary_arrow_table))

        # add the serialized data first (build flatbuffer from back to front)
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
    def read_data_partitions_as_arrow_table(cls, path_to_directory, file_ext='arrow'):
        cls.logger.info('>>> reading binary data partitions into an arrow table')

        # initialize these variables, just in case there are no files in path_to_directory
        table_partitions = []

        for partition_file in os.listdir(path_to_directory):
            if not partition_file.endswith(f'.{file_ext}'): continue

            with open(os.path.join(path_to_directory, partition_file), 'rb') as input_handle:
                table_partitions.append(arrow_table_from_binary(input_handle.read()))

        cls.logger.info(f'<<< read {len(table_partitions)} partitions')
        return pyarrow.concat_tables(table_partitions)

    @classmethod
    def read_data_file_as_arrow_table(cls, path_to_infile):
        cls.logger.info('>>> reading binary data file into an arrow table')

        with open(path_to_infile, 'rb') as input_handle:
            data_table = arrow_table_from_binary(input_handle.read())

        cls.logger.info('<<< data read into arrow table')

        return data_table

    @classmethod
    def read_data_file_as_flatbuffer(cls, path_to_infile):
        cls.logger.info('>>> reading binary data file into a flatbuffer')
        with open(path_to_infile, 'rb') as input_handle:
            binary_data = input_handle.read()

        flatbuffer_obj = SkyhookFlatbufferMeta.from_binary_flatbuffer(binary_data)
        cls.logger.info('<<< data read into arrow table')

        return flatbuffer_obj

    @classmethod
    def read_skyhook_file(cls, path_to_infile):
        cls.logger.info('>>> reading skyhook file into a flatbuffer')
        with open(path_to_infile, 'rb') as input_handle:
            binary_data = input_handle.read()

        flatbuffer_obj = SkyhookFlatbufferMeta.from_binary_flatbuffer(binary_data[4:])
        cls.logger.info('<<< data read into arrow table')

        return int.from_bytes(binary_data[:4], byteorder='little'), flatbuffer_obj

    @classmethod
    def read_data_file_as_binary(cls, path_to_infile):
        cls.logger.info('>>> reading binary data file')
        with open(path_to_infile, 'rb') as input_handle:
            binary_data = input_handle.read()

        cls.logger.info('<<< data read as binary')
        return binary_data


class SkyhookFileWriter(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def write_to_arrow(cls, data_wrapper, path_to_outfile):
        data_schema = data_wrapper.table_schema().with_metadata(
            data_wrapper.schema_skyhook_metadata().to_byte_coercible()
        )

        cls.logger.info('>>> writing data in single arrow file')
        with open(path_to_outfile, 'wb') as arrow_handle:
            batch_writer = pyarrow.RecordBatchFileWriter(arrow_handle, data_schema)

            for record_batch in data_wrapper.as_arrow_table(data_schema).to_batches():
                batch_writer.write_batch(record_batch)

        cls.logger.info('<<< data written')

    @classmethod
    def write_to_parquet(cls, data_wrapper, path_to_outfile):
        data_schema = data_wrapper.table_schema().with_metadata(
            data_wrapper.schema_skyhook_metadata().to_byte_coercible()
        )

        cls.logger.info('>>> writing data in single parquet file')
        pyarrow.parquet.write_table(data_wrapper.as_arrow_table(data_schema), path_to_outfile)

        cls.logger.info('<<< data written')

    @classmethod
    def write_partitions_to_flatbuffer(cls, data_wrapper, output_dir,
                                       batch_size=100, file_ext='skyhook'):
        # tbl_part is table partition
        for ndx, tbl_part in data_wrapper.batched_table_partitions(batch_size):
            # Generate path to binary file based on table partition metadata
            path_to_blob = os.path.join(
                output_dir,
                '{:05d}-{}.{}'.format(
                    ndx,
                    tbl_part.column_names[0],
                    file_ext
                )
            )

            # Serialize flatbuffer structure to disk
            fb_meta_wrapped_binary = SkyhookFlatbufferMeta.binary_from_arrow_binary(
                arrow_binary_from_table(tbl_part)
            )
            with open(path_to_blob, 'wb') as blob_handle:
                blob_handle.write(fb_meta_wrapped_binary)

    @classmethod
    def write_partitions_to_arrow(cls, data_wrapper, output_dir,
                                  batch_size=100, file_ext='arrow'):
        cls.logger.info('>>> writing data partitions into directory of arrow files')

        # tbl_part is table partition
        for ndx, tbl_part in data_wrapper.batched_table_partitions(batch_size=batch_size):
            # Generate path to binary file based on table partition metadata
            path_to_blob = os.path.join(
                output_dir,
                '{:05d}-{}.{}'.format(
                    ndx,
                    tbl_part.column_names[0],
                    file_ext
                )
            )

            # Serialize flatbuffer structure to disk
            arrow_binary_data = arrow_binary_from_table(tbl_part)
            with open(path_to_blob, 'wb') as blob_handle:
                blob_handle.write(arrow_binary_data)

        cls.logger.info('<<< data written')

    @classmethod
    def write_partitions_to_parquet(cls, data_wrapper, path_to_directory,
                                    batch_size=1000, file_ext='parquet'):
        cls.logger.info('>>> writing data partitions into directory of parquet files')

        # tbl_part is table partition
        for ndx, tbl_part in data_wrapper.batched_table_partitions(batch_size=batch_size):
            path_to_parquet_file = os.path.join(
                path_to_directory,
                '{:05d}-{}.{}'.format(
                    ndx,
                    tbl_part.column_names[0],
                    file_ext
                )
            )

            pyarrow.parquet.write_table(tbl_part, path_to_parquet_file)

        cls.logger.info('<<< data written')
