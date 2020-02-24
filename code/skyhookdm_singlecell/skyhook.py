from collections import namedtuple, OrderedDict


# ------------------------------
# Module-level Variables

skyhook_metadata_attributes = [
    'METADATA_SKYHOOK_VERSION'       ,
    'METADATA_DATA_SCHEMA_VERSION'   ,
    'METADATA_DATA_STRUCTURE_VERSION',
    'METADATA_DATA_FORMAT_TYPE'      ,
    'METADATA_DATA_SCHEMA'           ,
    'METADATA_DB_SCHEMA'             ,
    'METADATA_TABLE_NAME'            ,
    'METADATA_NUM_ROWS'              ,
]

column_schema_attributes = [
    'col_id'     ,
    'type'       ,
    'is_key'     ,
    'is_nullable',
    'col_name'   ,
]


# ------------------------------
# Classes
class SkyhookMetadata(namedtuple('SkyhookMetadata', skyhook_metadata_attributes)):
    def to_bytes(self):
        def bytes_from_val(val):
            if type(val) is int:
                return val.to_bytes(4, byteorder='big')

            elif type(val) is str:
                return bytes(val.encode('utf-8'))

            raise TypeError('Unable to convert value to bytes: {}({})'.format(val, type(val)))


        dict_entries = []
        for dict_key, dict_val in self._asdict().items():
            val_bytes = bytes_from_val(dict_val)

            dict_entries.append((dict_key, dict_val))

        return OrderedDict(dict_entries)


class ColumnSchema(namedtuple('ColumnSchema', column_schema_attributes)):
    def __str__(self):
        return ' '.join([str(field_val) for field_name, field_val in self._asdict().items()])


# Enums
class EnumMetaClass(type):

    def __new__(cls, name, bases, classdict_base):
        classdict_initialized = dict(classdict_base)

        # set the values of each slot, in the way that C/C++ sets enum values
        for enum_ndx, enum_name in enumerate(classdict_base['__enum_names__'], start=1):
            classdict_initialized[enum_name] = enum_ndx

        # always set first and last, so we generalize
        classdict_initialized['SDT_FIRST'] = 1
        classdict_initialized['SDT_LAST']  = enum_ndx

        # return the results of calling the superclass's __new__, but with the modified dictionary
        return type.__new__(cls, name, bases, classdict_initialized)


class DataTypes(object, metaclass=EnumMetaClass):
    """
    A class that represents Skyhook's SDT enum.
    """

    __slots__      = ()
    __enum_names__ = [
        'SDT_INT8' , 'SDT_INT16' , 'SDT_INT32' , 'SDT_INT64' ,
        'SDT_UINT8', 'SDT_UINT16', 'SDT_UINT32', 'SDT_UINT64',
        'SDT_CHAR' , 'SDT_UCHAR' , 'SDT_BOOL'  ,
        'SDT_FLOAT', 'SDT_DOUBLE',
        'SDT_DATE' , 'SDT_STRING',
    ]


class FormatTypes(object, metaclass=EnumMetaClass):
    """
    A class that represents Skyhook's SFT enum.
    """

    __slots__      = ()
    __enum_names__ = [
        'SFT_FLATBUF_FLEX_ROW'   , 'SFT_FLATBUF_UNION_ROW',
        'SFT_FLATBUF_UNION_COL'  , 'SFT_FLATBUF_CSV_ROW'  ,
        'SFT_ARROW'              , 'SFT_PARQUET'          ,
        'SFT_PG_TUPLE', 'SFT_CSV', 'SFT_PG_BINARY'        ,
        'SFT_HDF5'               , 'SFT_JSON'             ,
    ]


# Simple Types
class SentinelType(object):
    __slots__ = ()


class KeyColumn(SentinelType):
    NOT_KEY = 0
    KEY     = 1


class NullableColumn(SentinelType):
    NOT_NULLABLE = 0
    NULLABLE     = 1
