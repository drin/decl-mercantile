import logging

import flatbuffers
import pyarrow

# classes
from skyhookdm_singlecell.Tables.FB_Meta import FB_Meta

# functions
from skyhookdm_singlecell.Tables.FB_Meta import (FB_MetaStart, FB_MetaEnd,
                                                 FB_MetaAddBlobFormat,
                                                 FB_MetaAddBlobData,
                                                 FB_MetaAddBlobSize,
                                                 FB_MetaAddBlobDeleted,
                                                 FB_MetaAddBlobOrigOff,
                                                 FB_MetaAddBlobOrigLen,
                                                 FB_MetaAddBlobCompression)

class SkyhookFlatbuffer(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_binary(cls, flatbuffer_binary, offset=0):
        return cls(FB_Meta.GetRootAsFB_Meta(flatbuffer_binary, offset))

    def __init__(self, flatbuffer_obj, **kwargs):
        super().__init__(**kwargs)

        self.fb_obj = flatbuffer_obj

    def is_arrow_type(self):
        return self.fb_obj.BlobFormat() == 1

    def get_size(self):
        return self.fb_obj.BlobDataLength()

    def get_data_as_arrow(self):
        data_blob          = self.fb_obj.BlobDataAsNumpy()

        stream_reader      = pyarrow.ipc.open_stream(data_blob.tobytes())
        deserialized_table = pyarrow.Table.from_batches(
            list(stream_reader),
            schema=stream_reader.schema
        )

        return deserialized_table

class SkyhookFlatbufferGeneExpression(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_arrow_buffer(cls, arrow_buffer):
        return cls(initial_byte_count=arrow_buffer.size).wrap_arrow_buffer(arrow_buffer)

    def __init__(self, initial_byte_count=0, **kwargs):
        super().__init__(**kwargs)

        self.initial_byte_count = initial_byte_count
        self.builder            = flatbuffers.Builder(initial_byte_count)

    def wrap_arrow_buffer(self, arrow_buffer):
        """
        Write object, in-order
        """

        # we want to serialize the data first (build vector from back to front)
        wrapped_data_blob = self.builder.CreateByteVector(arrow_buffer.to_pybytes())

        FB_MetaStart(self.builder)

        FB_MetaAddBlobFormat(self.builder, 1)
        FB_MetaAddBlobData(self.builder, wrapped_data_blob)
        FB_MetaAddBlobSize(self.builder, arrow_buffer.size)
        FB_MetaAddBlobDeleted(self.builder, False)
        FB_MetaAddBlobOrigOff(self.builder, 0)
        FB_MetaAddBlobOrigLen(self.builder, arrow_buffer.size)
        FB_MetaAddBlobCompression(self.builder, 0)

        self.builder.Finish(FB_MetaEnd(self.builder))

        return self.builder.Output()
