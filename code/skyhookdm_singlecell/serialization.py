import logging

import flatbuffers

from skyhookdm_singlecell.Tables import FB_Meta

class SkyhookFlatbufferGeneExpression(object):
    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    @classmethod
    def from_arrow_buffer(cls, arrow_buffer):
        return cls(arrow_buffer.size).wrap_arrow_buffer(arrow_buffer)

    def __init__(self, initial_byte_count=0, **kwargs):
        super().__init__(**kwargs)

        self.initial_byte_count = initial_byte_count
        self.builder            = flatbuffers.Builder(initial_byte_count)

    def serialize_records(self, arrow_buffer):
        """
        Write element vector, in reverse order
        """

        FB_Meta.FB_MetaStartBlobDataVector(self.builder, arrow_buffer.size)

        for arrow_buffer_byte in reversed(arrow_buffer):
            self.builder.PrependByte(arrow_buffer_byte)

        return builder.EndVector(arrow_buffer.size)

    def wrap_arrow_buffer(self, arrow_buffer):
        """
        Write object, in-order
        """

        # we want to serialize the data first (build vector from back to front)
        wrapped_data_blob = self.serialize_records(arrow_buffer)

        FB_Meta.FB_MetaStart(self.builder)

        FB_Meta.FB_MetaAddBlobFormat(self.builder, 1)
        FB_Meta.FB_MetaAddBlobData(self.builder, wrapped_data_blob)
        FB_Meta.FB_MetaAddBlobSize(self.builder, arrow_buffer.size)
        FB_Meta.FB_MetaAddBlobDeleted(self.builder, False)
        FB_Meta.FB_MetaAddBlobOrigOff(self.builder, 0)
        FB_Meta.FB_MetaAddBlobOrigLen(self.builder, arrow_buffer.size)
        FB_Meta.FB_MetaAddBlobCompression(self.builder, 0)

        self.builder.Finish(FB_Meta.FB_MetaEnd(builder))

        return self.builder.Output()

    def deserialize(self, flatbuffer_obj, blob_byte_offset=0):
        skyhook_blob_wrapper = FB_Meta.GetRootAsFB_Meta(flatbuffer_obj, blob_byte_offset)

        skyhook_blob_wrapper.BlobDataLength()
        skyhook_blob_wrapper.BlobFormat()
