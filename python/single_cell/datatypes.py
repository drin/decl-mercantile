import os
import sys

import pyarrow

from collections import namedtuple

from single_cell.parsers import MatrixParser, GeneListParser, CellIDParser

debug = True

# ------------------------------
# Constants

# To represent Skyhook enum: SDT
skyhook_data_types = {
    data_type_name: enum_ndx
    for enum_ndx, data_type_name in enumerate(
         [
             'SDT_INT8' , 'SDT_INT16' , 'SDT_INT32' , 'SDT_INT64' ,
             'SDT_UINT8', 'SDT_UINT16', 'SDT_UINT32', 'SDT_UINT64',
             'SDT_CHAR' , 'SDT_UCHAR' , 'SDT_BOOL'  ,
             'SDT_FLOAT', 'SDT_DOUBLE',
             'SDT_DATE' , 'SDT_STRING',
         ]
        ,start=1
    )
}

skyhook_data_types['SDT_FIRST'] = skyhook_data_types['SDT_INT8']
skyhook_data_types['SDT_LAST']  = skyhook_data_types['SDT_STRING']


# To represent Skyhook enum: SFT
skyhook_format_types = {
    data_type_name: enum_ndx
    for enum_ndx, data_type_name in enumerate(
         [
             'SFT_FLATBUF_FLEX_ROW'   , 'SFT_FLATBUF_UNION_ROW',
             'SFT_FLATBUF_UNION_COL'  , 'SFT_FLATBUF_CSV_ROW'  ,
             'SFT_ARROW'              , 'SFT_PARQUET'          ,
             'SFT_PG_TUPLE', 'SFT_CSV', 'SFT_PG_BINARY'        ,
             'SFT_HDF5'               , 'SFT_JSON'             ,
         ]
        ,start=1
    )
}


# ------------------------------
# Simple classes
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

def skyhook_dataschema_from_gene_expr(gene_expr_obj, delim='\n'):
    schema_format = 'col_id type is_key is_nullable col_name;'

    # 'col_id type is_key is_nullable col_name;'
    # for type, I assume that 12 is the pos of the enum 'SDT_FLOAT'
    # for now, everything is nullable, nothing is a key, and everything's a float
    return '\n'.join([
        '{} {} {} {} {}'.format(
            cell_ndx,
            skyhook_data_types['SDT_FLOAT'],
            0,
            1,
            cell_id
        )
        for cell_ndx, cell_id in enumerate(gene_expr_obj.cells)
    ])


# ------------------------------
# Functions
def schema_from_gene_expr(gene_expr_obj):
    """
    Function that defines a pyarrow schema using a GeneExpression object.
    """

    gene_expr_metadata = SkyhookMetadata(
        0,
        0,
        0,
        skyhook_format_types['SFT_ARROW'],
        skyhook_dataschema_from_gene_expr(gene_expr_obj),
        '',
        'gene_expression',
        gene_expr_obj.expression.shape[0]
    )

    return pyarrow.schema((
        (cell_id, pyarrow.float64())
        for cell_id in gene_expr_obj.cells
    ))

def gene_expr_as_recordbatch(gene_expr_schema, gene_expr_obj):
    expr_data            = gene_expr_obj.expression.toarray()
    row_count, col_count = gene_expr_obj.expression.shape

    return pyarrow.RecordBatch.from_arrays(
         [
             pyarrow.array(expr_data[:, barcode_ndx])
             for barcode_ndx in range(col_count)
         ]
        ,schema=gene_expr_schema
    )


# ------------------------------
# utility functions
def write_recordbatches(path_to_outfile, batch_size, batch_offset, data_schema, data_recordbatch):
    batch_ndx = 0

    print('batch size: {}'.format(batch_size))
    print('batch offset: {}'.format(batch_offset))

    with open('{}'.format(path_to_outfile), 'wb') as output_handle:
        batch_writer = pyarrow.RecordBatchStreamWriter(output_handle, schema=data_schema)

        while batch_ndx + batch_offset < data_recordbatch.num_rows:

            if debug:
                sys.stdout.write('\rwriting batch: {:0d}'.format(batch_ndx))
                sys.stdout.flush()

            batch_writer.write_batch(data_recordbatch.slice(offset=batch_ndx, length=batch_size))
            batch_ndx += batch_offset

        # have to serialize the remainder batch (num_rows % batch_offset)
        if data_recordbatch.num_rows > batch_ndx:
            sys.stdout.write('\rwriting batch: {:0d}'.format(batch_ndx))
            sys.stdout.flush()

            batch_writer.write_batch(data_recordbatch.slice(offset=batch_ndx, length=batch_size))

        print('--- batches written')


# ------------------------------
# Classes
class GeneExpression(object):

    @classmethod
    def from_mtx(cls, path_to_mtx_root):
        """
        Factory function for reading gene expression matrix. Uses scanpy's read() function to read
        the gene expression data in mtx format, and then uses pandas to parse genes and cell IDs
        into dataframes to be used for obs and var dimensions of the gene expression anndata
        object.

        :path_to_mtx_root: Path to a directory containing: a gene expression matrix in mtx format,
        named 'matrix.mtx'; a list of cell IDs, named 'cells.tsv'; and a list of
        gene symbols (aka gene names), named 'genes.tsv'. 
        """

        # TODO this should eventually check for plaintext files too, but for now it assumes gzipped
        path_to_matrix = os.path.join(path_to_mtx_root, 'matrix.mtx')

        if not os.path.isfile(path_to_matrix):
            err_msg = 'Expected matrix.mtx in mtx root: {}\n'.format(path_to_mtx_root)

            sys.stderr.write(err_msg)
            sys.exit(err_msg)


        expr_matrix =   MatrixParser.from_path(os.path.join(path_to_mtx_root, 'matrix.mtx'))
        gene_list   = GeneListParser.from_path(os.path.join(path_to_mtx_root, 'genes.tsv'))
        cells       =   CellIDParser.from_path(os.path.join(path_to_mtx_root, 'cells.tsv'))

        return cls(expr_matrix.tocsc(), gene_list, cells)

    def __init__(self, expr_matrix, gene_list, cells, **kwargs):
        super().__init__(**kwargs)

        self.expression = expr_matrix
        self.genes      = gene_list
        self.cells      = cells
