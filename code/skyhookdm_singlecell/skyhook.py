from collections import namedtuple, OrderedDict


# ------------------------------
# Module-level Variables

skyhook_metadata_attributes = [
    'skyhook_version'       ,
    'data_schema_version'   ,
    'data_structure_version',
    'data_format_type'      ,
    'data_schema'           ,
    'db_schema'             ,
    'table_name'            ,
    'num_rows'              ,
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
    def to_byte_coercible(self):
        def bytes_from_val(val):
            if type(val) is int:
                return val.to_bytes(4, byteorder='little')

            elif type(val) is str:
                return bytes(val.encode('utf-8'))

            raise TypeError('Unable to convert value to bytes: {}({})'.format(val, type(val)))

        return OrderedDict([
            (dict_key, bytes_from_val(dict_val))
            for dict_key, dict_val in self._asdict().items()
        ])


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

    def from_value(self, format_type_val):
        if not format_type_val: return None

        return self.__enum_names__[format_type_val - 1]


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
        'SFT_FLATBUF_FLEX_ROW' ,
        'SFT_FLATBUF_UNION_ROW',
        'SFT_FLATBUF_UNION_COL',
        'SFT_FLATBUF_CSV_ROW'  ,
        'SFT_ARROW'            ,
        'SFT_PARQUET'          ,
        'SFT_PG_TUPLE'         ,
        'SFT_CSV'              ,
        'SFT_JSON'             ,
        'SFT_PG_BINARY'        ,
        'SFT_PYARROW_BINARY'   ,
        'SFT_HDF5'             ,
        'SFT_EXAMPLE_FORMAT'   ,
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
