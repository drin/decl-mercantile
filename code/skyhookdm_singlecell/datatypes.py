# core packages
import os
import sys

# dependencies
import pyarrow

# sub-modules of this package
import skyhook

from skyhookdm_singlecell import (__skyhook_version__            ,
                                  __skyhook_data_schema_version__,
                                  __skyhook_data_struct_version__,)

from parsers import MatrixParser, GeneListParser, CellIDParser


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

        expr_matrix = MatrixParser.from_path(os.path.join(path_to_mtx_root, 'matrix.mtx'))
        gene_list   = GeneListParser.from_path(os.path.join(path_to_mtx_root, 'genes.tsv'))
        cells       = CellIDParser.from_path(os.path.join(path_to_mtx_root, 'cells.tsv'))

        return cls(expr_matrix.tocsc(), gene_list, cells)

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
