# core packages
import sys

# dependencies
import numpy
import scipy


def cpm_normalization(expression_data, target_sum=1e6):
    """
    Take the log(<count> + 1) for each matrix cell (for each barcode, for each gene) where each
    measured gene expression (matrix cell) is normalized by the proportion of total expressed
    genes (total count for the barcode or biological cell) to 1 million.

    In short, take the log(CPM + 1).
    """

    # normalization across expression in the single-cell by a million (default)
    rna_total_normalized = expression_data.sum(axis=0) / target_sum

    # normalize barcode counts by sum of normalized counts for that barcode
    rna_expression_normalized = expression_data / rna_total_normalized

    return (
        numpy.log1p(rna_expression_normalized)
             .astype(numpy.uint16, copy=False)
    )


# ------------------------------
# Classes
class Annotation(object):
    __slots__ = ('headers', 'annotations')

    def __init__(self, headers, annotations):
        self.headers     = headers
        self.annotations = annotations

    def __getitem__(self, col_name):
        return self.annotations[:, numpy.where(self.headers == col_name)]

    @property
    def shape(self):
        """Shape of expression matrix."""
        return self.annotations.shape

    def astype(self, dtype):
        self.annotations = self.annotations.astype(dtype=dtype, copy=False)

        return self

    def columns(self, from_col_ndx=None, to_col_ndx=None):
        return self.headers[from_col_ndx:to_col_ndx]

    def data(self, from_col=None, to_col=None, from_row=None, to_row=None):
        return self.annotations[from_row:to_row, from_col:to_col]

    def data_as_array(self, from_col=None, to_col=None, from_row=None, to_row=None):
        is_sparse_matrix = (
               isinstance(self.annotations, scipy.sparse.coo_matrix)
            or isinstance(self.annotations, scipy.sparse.csc_matrix)
        )

        if (is_sparse_matrix):
            data_array = self.annotations[from_row:to_row, from_col:to_col].toarray()

            if data_array.shape[1] == 1: return data_array[:,0]
            return data_array

        elif isinstance(self.annotations, numpy.matrix):
            return self.annotations[from_row:to_row, from_col:to_col].getA()

        elif isinstance(self.annotations, numpy.ndarray):
            return self.annotations[from_row:to_row, from_col:to_col]

        sys.exit(f'Not sure how to convert {type(self.annotations)} to {type(numpy.ndarray)}')


class GeneExpression(object):

    def __init__(self, expression_matrix, genes, cells, **kwargs):
        super().__init__(**kwargs)

        # Core data
        self.genes      = genes
        self.cells      = cells
        self.expression = expression_matrix

        if isinstance(expression_matrix, scipy.sparse.coo_matrix):
            self.expression = expression_matrix.tocsc()

        self.pipeline_functions = []

    @property
    def shape(self):
        """Shape of expression matrix."""
        return self.expression.shape

    def astype(self, dtype):
        self.expression = self.expression.astype(dtype=dtype, copy=False)

        return self

    def columns(self, from_col_ndx=None, to_col_ndx=None):
        return self.cells[from_col_ndx:to_col_ndx]

    def rows(self, from_row_ndx=None, to_row_ndx=None):
        return self.genes[from_row_ndx:to_row_ndx]

    def data(self, from_col=None, to_col=None, from_row=None, to_row=None):
        return self.expression[from_row:to_row, from_col:to_col]

    def data_as_array(self, from_col=None, to_col=None, from_row=None, to_row=None):
        data_array       = None
        is_sparse_matrix = (
               isinstance(self.expression, scipy.sparse.coo_matrix)
            or isinstance(self.expression, scipy.sparse.csc_matrix)
        )

        if (is_sparse_matrix):
            data_array = self.expression[from_row:to_row, from_col:to_col].toarray()

            if data_array.shape[1] == 1:
                data_array = data_array[:,0]

        elif isinstance(self.expression, numpy.matrix):
            data_array = self.expression[from_row:to_row, from_col:to_col].getA()

        elif isinstance(self.expression, numpy.ndarray):
            data_array = self.expression[from_row:to_row, from_col:to_col]

        else:
            sys.exit(f'Not sure how to convert {type(self.expression)} to numpy.ndarray')

        for fn_pipelined in self.pipeline_functions:
            data_array = fn_pipelined(data_array)

        return data_array

    def subsample_by_counts(self, gene_count=100, cell_count=10):
        if isinstance(self.expression, scipy.sparse.coo_matrix):
            self.expression = self.expression.tocsc().todense()

        self.expression = self.expression[:gene_count, :cell_count]
        self.genes      = self.genes[:gene_count]
        self.cells      = self.cells[:cell_count]

        return self

    def pipeline_add_cpm_normalization(self):
        self.pipeline_functions.append(cpm_normalization)

    def normalize_expression(self, target_sum=1e6):
        """
        By default, calls `cpm_normalization`, which computes the log(CPM + 1).
        """

        self.expression = cpm_normalization(self.expression, target_sum=target_sum)
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
