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
class Annotation(object):
    __slots__ = ('_headers', '_annotations')

    @classmethod
    def from_csv(cls, path_to_annotations, has_header=True, delim=','):
        """
        Returns annotations represented by a numpy.array object of shape:
            (# of parsed annotations, # of annotation columns)
        """

        with open(path_to_annotations, 'r') as ann_handle:
            # parse headers
            ann_headers = None
            if has_header:
                ann_headers = numpy.array(next(ann_handle).strip().split(delim), dtype=str)

            annotations = [
                numpy.array(ann_line.strip().split(delim), dtype=str)
                for ann_line in ann_handle
            ]

        # return tuple (header, data)
        return Annotation(ann_headers, numpy.array(annotations))


    def __init__(self, headers, annotations):
        self._headers     = headers
        self._annotations = annotations


    def __getitem__(self, col_name):
        return self._annotations[:, numpy.where(self._headers == col_name)]


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

    @classmethod
    def from_custom_mtx(cls, path_to_mtx_root):
        """
        Factory function for reading gene expression matrix that is very similar to
        GeneExpression.from_mtx() but applies a transformation to barcodes on ingest.

        :path_to_mtx_root: Path to a directory containing: a gene expression matrix in mtx format,
        named 'matrix.mtx'; a list of barcodes (aka Cell IDs), named 'barcodes.tsv'; and a list of
        gene symbols (aka gene names), named 'genes.tsv'. 
        """

        try:
            sample_name = os.path.basename(path_to_mtx_root).split('-')[1]
            closure_transform_barcode = cls.cell_id_from_barcode(sample_name)

        except Exception:
            err_msg = 'Unable to parse sample name from path: {}'.format(path_to_mtx_root)

            sys.stderr.write(err_msg)
            sys.exit(err_msg)

        gene_expr = MatrixParser.from_path(os.path.join(path_to_mtx_root, 'matrix.mtx'))
        gene_list = GeneListParser.from_path(os.path.join(path_to_mtx_root, 'genes.tsv'))
        barcodes  = CellIDParser.from_path(
            os.path.join(path_to_mtx_root, 'barcodes.tsv'),
            fn_transform=closure_transform_barcode
        )

        return cls(gene_expr, gene_list, barcodes)

    @classmethod
    def cell_id_from_barcode(cls, sample_name):
        """
        This is a utility function that returns a closure to be mapped onto a list of barcodes.
        """

        def closure_cell_id_from_barcode(cell_barcode):
            """
            To get 'Cell ID' we prefix the cell barcode with '<sample_name>_' and we remove the
            '-1' from the end.
            """

            return '{}_{}'.format(sample_name, cell_barcode.rstrip('-1'))

        return closure_cell_id_from_barcode

    def __init__(self, expr_matrix, gene_list, cells, **kwargs):
        super().__init__(**kwargs)

        self.expression = expr_matrix
        self.genes      = gene_list
        self.cells      = cells

    def normalize_expression(self, target_sum=1e6):
        """
        Note that log1p is called when the gene expression matrix is transposed to genes (obs;
        rows) x cells (var; cols) to match Bianca's code, but scanpy seems to assume that log1p is
        given a matrix of cells (obs; rows) x genes (var; cols). See the documentation on the first
        parameter here:
            https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.log1p.html#scanpy.pp.log1p
        """

        # if the gene expression data is in a sparse, coordinate matrix (mtx format), then we can
        # convert to a Compressed, Sparse Column (csc) matrix to be able to do operations on it.
        if isinstance(self._expr_data, scipy.sparse.coo_matrix):
            self._expr_data = self._expr_data.tocsc()

        # Take the log(<count> + 1) for each matrix cell (for each barcode, for each gene) where
        # each measured gene expression (matrix cell) is normalized by the proportion of total
        # expressed genes (total count for the barcode or biological cell) to 1 million.
        # In short, take the log(CPM + 1).
        normalized_counts = numpy.log1p(
              # normalize each count for a barcode by the total normalized count for that barcode
              self._expr_data
            / (
                  # normalize total count for each barcode by a million
                  self._expr_data.sum(axis=0) / target_sum
              )
        )

        return self


    def filter_cells_by_annotation(self, ann_data, count_threshold=50):
        """
        :ann_data: is expected to be an Annotation object.
        """

        if not type(ann_data) is Annotation:
            err_msg = 'Expected Annotation object, got "{}"'.format(type(ann_data))

            sys.stderr.write(err_msg)
            sys.exit(err_msg)

        # for now we just assume the annotation file contains a proper superset of the cells in the
        # gene expression matrix
        cell_annotations = Annotation(
             ann_data._headers
            ,ann_data[numpy.where(
                 numpy.isin(ann_data['cell'], self._barcodes)
             )]
        )

        filtered_cell_types  = [
            cell_type
            for cell_type in numpy.unique(cell_annotations['cell_ontology_class'])
            if numpy.sum(cell_annotations['cell_ontology_class'] == cell_type) > count_threshold
        ]

        # TODO: verify that barcodes_supported_by_cell_type properly co-indexes barcodes
        barcodes_supported_by_cell_type = numpy.where(
            numpy.isin(annotated_cell_types, filtered_cell_types)
        )

        # use filtered barcodes to filter matrix barcodes and columns of the expression
        self._expr_data = self._expr_data[:,barcodes_supported_by_cell_type]
        self._barcodes  =  self._barcodes[barcodes_supported_by_cell_type]

        return self


    def compute_expression_by_cell_type(self, ann_data, cell_type):
        """
        A function that will:
            * compute expression for each cell type (cell_ontology_class)
            * compute expression for N, N/2, N/3, N/4, N/5 statistical samples
              of the gene expression dataset
            * should return information necessary for tsne visualization
        """

        # get the annotations that have cell IDs in our barcode dataset
        cell_annotations = Annotation(
             ann_data._headers
            ,ann_data[numpy.where(
                 numpy.isin(ann_data['cell'], self._barcodes)
             )]
        )

        # get the annotation entries that match the desired cell_type
        col_indices = numpy.where(
            numpy.isin(cell_type, cell_annotations['cell_ontology_class'])
        )

        # get the working set of barcodes and the expression matrix
        workset_barcodes  = self._barcodes[col_indices]
        workset_expr_data = self._expr_data[:,col_indices]

        # TODO
        # tpm_avg = numpy.
