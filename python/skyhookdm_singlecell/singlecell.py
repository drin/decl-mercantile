# core packages
import sys

# dependencies
import numpy
import scipy


# ------------------------------
# Classes
class Annotation(object):
    __slots__ = ('headers', 'annotations')

    def __init__(self, headers, annotations):
        self.headers     = headers
        self.annotations = annotations

    def __getitem__(self, col_name):
        return self.annotations[:, numpy.where(self.headers == col_name)]


class GeneExpression(object):

    def __init__(self, expression_matrix, genes, cells, **kwargs):
        super().__init__(**kwargs)

        # Core data
        self.expression = expression_matrix
        self.genes      = genes
        self.cells      = cells

    def subsample_by_counts(self, gene_count=100, cell_count=10):
        if isinstance(self.expression, scipy.sparse.coo_matrix):
            self.expression = self.expression.tocsc().todense()

        self.expression = self.expression[:gene_count, :cell_count]
        self.genes      = self.genes[:gene_count]
        self.cells      = self.cells[:cell_count]

        return self

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
        if isinstance(self.expression, scipy.sparse.coo_matrix):
            self.expression = self.expression.tocsc()

        # Take the log(<count> + 1) for each matrix cell (for each barcode, for each gene) where
        # each measured gene expression (matrix cell) is normalized by the proportion of total
        # expressed genes (total count for the barcode or biological cell) to 1 million.
        # In short, take the log(CPM + 1).
        normalized_counts = numpy.log1p(
            # normalize each count for a barcode by the total normalized count for that barcode
            self.expression / (
                # normalize total count for each barcode by a million
                self.expression.sum(axis=0) / target_sum
            )
        )

        self.expression = normalized_counts

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
             ann_data.headers
            ,ann_data[numpy.where(
                 numpy.isin(ann_data['cell'], self.cells)
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
        self.expression = self.expression[:,barcodes_supported_by_cell_type]
        self.cells      = self.cells[barcodes_supported_by_cell_type]

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
             ann_data.headers
            ,ann_data[numpy.where(
                 numpy.isin(ann_data['cell'], self.cells)
             )]
        )

        # get the annotation entries that match the desired cell_type
        col_indices = numpy.where(
            numpy.isin(cell_type, cell_annotations['cell_ontology_class'])
        )

        # get the working set of barcodes and the expression matrix
        workset_barcodes  = self.cells[col_indices]
        workset_expr_data = self.expression[:,col_indices]

        workset_barcodes, workset_expr_data

        # TODO
        # tpm_avg = numpy.
