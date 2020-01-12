from collections import namedtuple


# ------------------------------
# Module-level Variables
debug = True

SkyhookMetadata = namedtuple('SkyhookMetadata', [
    'METADATA_SKYHOOK_VERSION'       ,
    'METADATA_DATA_SCHEMA_VERSION'   ,
    'METADATA_DATA_STRUCTURE_VERSION',
    'METADATA_DATA_FORMAT_TYPE'      ,
    'METADATA_DATA_SCHEMA'           ,
    'METADATA_DB_SCHEMA'             ,
    'METADATA_TABLE_NAME'            ,
    'METADATA_NUM_ROWS'              ,
])


# ------------------------------
# Classes
class EnumMetaClass(type):

    def __new__(cls, name, bases, classdict_base):
        classdict_initialized = dict(classdict_base)

        # set the values of each slot, in the way that C/C++ sets enum values
        for enum_ndx, slot_name in enumerate(classdict_base['__slots__'], start=1):
            classdict_initialized[slot_name] = enum_ndx

        # always set first and last, so we generalize
        classdict_initialized['SDT_FIRST'] = 1
        classdict_initialized['SDT_LAST']  = enum_ndx

        # return the results of calling the superclass's __new__, but with the modified dictionary
        return type.__new__(cls, name, bases, classdict_initialized)


class DataTypes(object, metaclass=EnumMetaClass):
    """
    A class that represents Skyhook's SDT enum.
    """

    __slots__ = [
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

    __slots__ = [
        'SFT_FLATBUF_FLEX_ROW'   , 'SFT_FLATBUF_UNION_ROW',
        'SFT_FLATBUF_UNION_COL'  , 'SFT_FLATBUF_CSV_ROW'  ,
        'SFT_ARROW'              , 'SFT_PARQUET'          ,
        'SFT_PG_TUPLE', 'SFT_CSV', 'SFT_PG_BINARY'        ,
        'SFT_HDF5'               , 'SFT_JSON'             ,
    ]


# ------------------------------
# Functions
def skyhook_dataschema_from_gene_expr(gene_expr_obj, delim='\n'):
    # TODO: Not yet used
    # schema_format = 'col_id type is_key is_nullable col_name;'

    # 'col_id type is_key is_nullable col_name;'
    # for type, I assume that 12 is the pos of the enum 'SDT_FLOAT'
    # for now, everything is nullable, nothing is a key, and everything's a float
    return '\n'.join([
        '{} {} {} {} {}'.format(
            cell_ndx,
            DataTypes.SDT_FLOAT,
            0,
            1,
            cell_id
        )
        for cell_ndx, cell_id in enumerate(gene_expr_obj.cells)
    ])
