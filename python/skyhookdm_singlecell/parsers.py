import os
import sys
import gzip

import numpy
import scipy.io

# classes
from skyhookdm_singlecell.datatypes import GeneExpression, Annotation

# functions
from skyhookdm_singlecell.util import normalize_str


class GeneListParser(object):

    __slots__ = ()

    @classmethod
    def from_handle(cls, genelist_handle, has_header=True, delim='\t'):
        """
        This function iterates over each row in :genelist_handle:, building a numpy array
        containing each field of the row. Then, each row is accumulated into a python list to be
        returned as a numpy array (records from genelist file) or numpy arrays (fields for each
        genelist row).

        :genelist_handle: (file-like)              A path to the metadata for genes.
        :has_header:      (bool)                   A flag, determining whether to skip the first
                                                   line or not.
        :delim:           (string) (default: '\t') The delimiter to use to distinguish between
                                                   attributes (fields) of each gene (for each row
                                                   in the genelist file).
        """

        if not genelist_handle or genelist_handle.closed:
            err_msg = 'Gene list file handle is closed or never opened.'

            sys.stderr.write(f'{err_msg}\n')
            sys.exit(err_msg)

        if has_header: next(genelist_handle)

        gene_list = [
            numpy.array(normalize_str(gene_line).strip().strip('"').split(delim), dtype=str)
            for gene_line in genelist_handle
        ]

        return numpy.array(gene_list, dtype=str)


class CellIDParser(object):

    __slots__ = ()

    @classmethod
    def from_handle(cls, cell_list_handle, has_header=True, delim='\t', fn_transform=None):
        """
        This function currently hardcodes the column to find the primary key for cells. For the
        MOCA dataset, the file was single column and straightforward. HCA metadata is incredibly
        more complex, but goes through the effort of uniquely identifying a cell using a "cellkey"
        attribute. I do not know how uniqueness of cellkey is determined.

        Currently, no validation of the file contents is done.
        """

        if not cell_list_handle or cell_list_handle.closed:
            err_msg = 'Cell list file handle is closed or never opened.'

            sys.stderr.write(f'{err_msg}\n')
            sys.exit(err_msg)

        cellkey_col_ndx = 0

        with cell_list_handle:
            if has_header:
                next(cell_list_handle)
                '''
                Currently, not used, so just skip the header if it exists
                metadata_columns = numpy.array(
                    next(cell_list_handle).strip() .split(delim),
                    dtype=str
                )
                '''

            # eventually migrate to a dynamic access such as this
            # target_columns = numpy.where(metadata_columns in ['barcode',])

            cell_ids = [
                normalize_str(cell_record).strip().strip('"').split(delim)[cellkey_col_ndx]
                for cell_record in cell_list_handle
            ]

        # return gene list as a numpy array of strings
        if fn_transform:
            return numpy.array(list(map(fn_transform, cell_ids)), dtype=str)

        return numpy.array(cell_ids, dtype=str)


# ------------------------------
# file format parsers
class AnnotationCSVParser(object):

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


class GeneExpressionMatrixParser(object):

    __slots__ = ()

    # ------------------------------
    # Convenience functions used by the parsing API
    @classmethod
    def open_component_from_dir(cls, path_to_mtx_root, component_name):
        # make some attempts to find the matrix file

        uncompressed_path = os.path.join(path_to_mtx_root, component_name)
        compressed_path   = os.path.join(path_to_mtx_root, f'{component_name}.gz')

        if os.path.isfile(uncompressed_path):
            return open(uncompressed_path, 'r')

        elif os.path.isfile(compressed_path):
            return gzip.open(compressed_path, 'rb')

        # error checking for matrix file
        else:
            err_msg = f'Expected {component_name}[.gz] in mtx root: {path_to_mtx_root}\n'

            sys.stderr.write(err_msg)
            sys.exit(err_msg)

        return None

    @classmethod
    def open_matrix_from_dir(cls, path_to_mtx_root):
        return cls.open_component_from_dir(path_to_mtx_root, 'matrix.mtx')

    @classmethod
    def open_genes_from_dir(cls, path_to_mtx_root):
        return cls.open_component_from_dir(path_to_mtx_root, 'genes.tsv')

    @classmethod
    def open_cells_from_dir(cls, path_to_mtx_root):
        return cls.open_component_from_dir(path_to_mtx_root, 'cells.tsv')

    @classmethod
    def open_features_from_dir(cls, path_to_mtx_root):
        """
        A component dataset that seems specific to recent versions of HCA datasets.
        """

        return cls.open_component_from_dir(path_to_mtx_root, 'features.tsv')

    @classmethod
    def open_barcodes_from_dir(cls, path_to_mtx_root):
        """
        A component dataset that is specially named in the MOCA project.
        """

        return cls.open_component_from_dir(path_to_mtx_root, 'barcodes.tsv')

    # ------------------------------
    # Auxiliary Functions
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

    # ------------------------------
    # Parsing API
    @classmethod
    def gene_expr_from_dir(cls, path_to_mtx_dir, has_header=True):
        """
        :path_to_mtx_dir: Path to directory containing matrix data, in MTX format. Files expected
                          in directory:
                              - 'matrix.mtx'
                              - 'genes.tsv'
                              - 'cells.tsv'
                              - 'features.tsv' (optional)

                          Files may be uncompressed, or gzipped ('.gz' suffix).
        """

        # Use scipy's function to read sparse matrix in MTX format.
        matrix_data = scipy.io.mmread(cls.open_matrix_from_dir(path_to_mtx_dir))

        gene_list = GeneListParser.from_handle(
            cls.open_genes_from_dir(path_to_mtx_dir),
            has_header=has_header
        )

        cell_list = CellIDParser.from_handle(
            cls.open_cells_from_dir(path_to_mtx_dir),
            has_header=has_header
        )

        return GeneExpression(
            expression_matrix=matrix_data.tocsc(),
            genes=gene_list,
            cells=cell_list
        )

    @classmethod
    def gene_expr_from_dir_custom(cls, path_to_mtx_dir, has_header=True):
        """
        :path_to_mtx_dir: Path to directory containing matrix data, in MTX format. Files expected
                          in directory:
                              - 'matrix.mtx'
                              - 'genes.tsv'
                              - 'cells.tsv'
                              - 'features.tsv' (optional)

                          Files may be uncompressed, or gzipped ('.gz' suffix).
        """

        try:
            sample_name = os.path.basename(path_to_mtx_dir).split('-')[1]
            closure_transform_barcode = cls.cell_id_from_barcode(sample_name)

        except Exception:
            err_msg = 'Unable to parse sample name from path: {}'.format(path_to_mtx_dir)

            sys.stderr.write(err_msg)
            sys.exit(err_msg)

        # Use scipy's function to read sparse matrix in MTX format.
        matrix_data = scipy.io.mmread(cls.open_matrix_from_dir(path_to_mtx_dir))

        gene_list = GeneListParser.from_handle(
            cls.open_genes_from_dir(path_to_mtx_dir),
            has_header=has_header
        )

        barcode_list = CellIDParser.from_handle(
            cls.open_barcodes_from_dir(path_to_mtx_dir),
            has_header=has_header,
            fn_transform=closure_transform_barcode
        )

        return GeneExpression(
            expression_matrix=matrix_data.tocsc(),
            genes=gene_list,
            cells=barcode_list
        )


class HCAMetadataParser(object):

    __slots__ = ()

    @classmethod
    def from_path(cls, path_to_metadata, has_header=True, delim='\t'):
        """
        This function assumes a many column file where the header line contains column names and no
        leading character (such as '#'). Currently, no validation of the file contents is done.
        """

        with open(path_to_metadata, 'r') as metadata_handle:
            # peel header, if necessary
            if has_header: next(metadata_handle)

            # TODO: this is for more dynamic parsing in the future, but isn't necessary now
            # if has_header:
            #     metadata_columns = numpy.array(
            #          next(metadata_handle).strip()
            #                               .split(delim)
            #         ,dtype=str
            #     )

            # target_columns = numpy.where(
            #     metadata_columns in [
            #         'project.project_core.project_short_name',
            #         'bundle_uuid',
            #         'file_uuid',
            #     ]
            # )

            # metadata = [
            #     metadata_line.strip().split(delim)
            #     for metadata_line in metadata_handle
            # ]

            # NOTE: we only need to peel one row to get these values because they're the same
            # throughout the metadata file
            metadata_fields = next(metadata_handle).split(delim)
            bundle_uuid     = metadata_fields[32]
            proj_name       = metadata_fields[-48]

        # return gene list as a numpy array of strings
        # return metadata_columns, numpy.array(metadata, dtype=str)
        return bundle_uuid, proj_name
