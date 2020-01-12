# core packages
import os
import sys

# dependencies
import numpy
import scipy
import pyarrow

# sub-modules of this package
from skyhookdm_singlecell import (__skyhook_version__            ,
                                  __skyhook_data_schema_version__,
                                  __skyhook_data_struct_version__,)

from skyhookdm_singlecell import skyhook

from skyhookdm_singlecell.parsers import MatrixParser, GeneListParser, CellIDParser


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

        expr_matrix = MatrixParser.from_path(os.path.join(path_to_mtx_root, 'matrix.mtx'))
        gene_list   = GeneListParser.from_path(os.path.join(path_to_mtx_root, 'genes.tsv'))
        cells       = CellIDParser.from_path(os.path.join(path_to_mtx_root, 'cells.tsv'))

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

        # Core data
        self.expression = expr_matrix
        self.genes      = gene_list
        self.cells      = cells

        # Skyhook metadata
        self.skyhook_schema = None
        self.db_schema      = ''
        self.table_name     = 'gene_expression'

    def clear_schema(self):
        self.skyhook_schema = None

    def to_arrow_recordbatch(self):
        expr_data            = self.expression.toarray()
        row_count, col_count = self.expression.shape

        return pyarrow.RecordBatch.from_arrays(
             [
                 pyarrow.array(expr_data[:, barcode_ndx])
                 for barcode_ndx in range(col_count)
             ]
            ,schema=self.skyhook_schema_arrow()
        )

    def skyhook_schema_arrow(self):
        if self.skyhook_schema is None:
            # create an arrow schema for every column (expression for a cell)
            # NOTE: pyarrow data type should match the data type used in skyhook data schema
            arrow_schema = pyarrow.schema((
                (cell_id, pyarrow.uint16())
                for cell_id in self.cells
            ))

            self.skyhook_schema = arrow_schema.with_metadata(self.skyhook_metadata()._asdict())

        return self.skyhook_schema

    def skyhook_metadata(self):
        """
        Return formatted Skyhook metadata that will be stored with the arrow schema information.
        """

        # the data schema is a list of column schemas
        skyhook_data_schema = '\n'.join([
            skyhook.ColumnSchema(
                cell_ndx                       ,
                skyhook.DataTypes.SDT_UINT16   ,
                skyhook.KeyColumn.NOT_KEY      ,
                skyhook.NullableColumn.NULLABLE,
                cell_id
            )
            for cell_ndx, cell_id in enumerate(self.cells)
        ])

        # return the Skyhook metadata tuple
        return skyhook.SkyhookMetadata(
            __skyhook_version__            ,
            __skyhook_data_schema_version__,
            __skyhook_data_struct_version__,
            skyhook.FormatTypes.SFT_ARROW  ,
            skyhook_data_schema            ,
            self.db_schema                 ,
            self.table_name                ,
            self.expression.shape[0]
        )

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

        # just for flake8 satisfaction; for now
        normalized_counts

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
            numpy.isin(cell_annotations['cell_ontology_class'], filtered_cell_types)
        )

        # use filtered barcodes to filter matrix barcodes and columns of the expression
        self._expr_data = self._expr_data[:,barcodes_supported_by_cell_type]
        self._barcodes  = self._barcodes[barcodes_supported_by_cell_type]

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

        workset_barcodes, workset_expr_data

        # TODO
        # tpm_avg = numpy.
