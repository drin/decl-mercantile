import os
import sys

import re
import gzip

import numpy
import scipy.io

# classes
from skyhookdm_singlecell.singlecell import GeneExpression, Annotation

# functions
from skyhookdm_singlecell.util import normalize_str


def normalize_file_record(record_as_str, delim='\t'):
    return (
        normalize_str(record_as_str).strip()
                                    .strip('"')
                                    .split(delim)
    )


# ------------------------------
# Query Parsers

class QueryParser(object):

    query_structure_template = r'SELECT\s+(.*)\s+FROM\s+(.*)'

    group_ndx_select = 1
    group_ndx_from   = 2

    @classmethod
    def parse_clause_select(cls, select_clause):
        projection_attrs = [
            projection_attr.strip()
            for projection_attr in select_clause.strip().split(',')
        ]

        if any(projection_attr == '*' for projection_attr in projection_attrs):
            return []

        return projection_attrs

    @classmethod
    def parse_clause_from(cls, from_clause):
        return [
            relation.strip()
            for relation in from_clause.strip().split(',')
        ]

    @classmethod
    def parse_clause_where(cls, where_clause):
        pass

    @classmethod
    def parse(cls, query_string):
        matched_query = re.search(cls.query_structure_template, query_string)

        projection_attrs = cls.parse_clause_select(matched_query.group(cls.group_ndx_select))
        query_relations  = cls.parse_clause_from(matched_query.group(cls.group_ndx_from))

        return {
            'projection': projection_attrs,
            'relations' : query_relations,
            'predicates': '',
            'aggregates': '',
        }


# ------------------------------
# Data Parsers

class GeneMetadataParser(object):

    __slots__ = ()

    @classmethod
    def from_handle(cls, gene_metadata_handle, has_header=True, delim='\t'):
        """
        This function iterates over each row in :gene_metadata_handle:, building a numpy array
        containing each field of the row. Then, each row is accumulated into a python list to be
        returned as a numpy array (records from genelist file) or numpy arrays (fields for each
        genelist row).

        :gene_metadata_handle: (file-like)              A path to the metadata for genes.
        :has_header:      (bool)                   A flag, determining whether to skip the first
                                                   line or not.
        :delim:           (string) (default: '\t') The delimiter to use to distinguish between
                                                   attributes (fields) of each gene (for each row
                                                   in the genelist file).
        """

        if not gene_metadata_handle or gene_metadata_handle.closed:
            err_msg = 'Gene list file handle is closed or was never opened.'

            sys.stderr.write(f'{err_msg}\n')
            sys.exit(err_msg)

        gene_header = ['featurekey']
        if has_header:
            gene_header = normalize_file_record(next(gene_metadata_handle), delim)

        gene_metadata = [
            numpy.array(normalize_file_record(gene_line, delim), dtype=str)
            for gene_line in gene_metadata_handle
        ]

        # ------------------------------
        # Returns a tuple: <array of headers>, <array of gene metadata rows>
        # Note that a row of gene metadata is, itself, an array
        return (
            numpy.array(gene_header  , dtype=str),
            numpy.array(gene_metadata, dtype=str)
        )


class CellMetadataParser(object):

    __slots__ = ()

    @classmethod
    def from_handle(cls, cell_metadata_handle, has_header=True, delim='\t', fn_transform=None):
        """
        This function currently hardcodes the column to find the primary key for cells. For the
        MOCA dataset, the file was single column and straightforward. HCA metadata is incredibly
        more complex, but goes through the effort of uniquely identifying a cell using a "cellkey"
        attribute. I do not know how uniqueness of cellkey is determined.

        Currently, no validation of the file contents is done.
        """

        if not cell_metadata_handle or cell_metadata_handle.closed:
            err_msg = 'Cell list file handle is closed or was never opened.'

            sys.stderr.write(f'{err_msg}\n')
            sys.exit(err_msg)

        cell_header = ['cellkey']
        if has_header:
            cell_header = normalize_file_record(next(cell_metadata_handle), delim)

        cell_metadata = [
            numpy.array(normalize_file_record(cell_record, delim), dtype=str)
            for cell_record in cell_metadata_handle
        ]

        # ------------------------------
        # Returns a tuple: <array of headers>, <array of cell metadata rows>
        # Note that a row of cell metadata is, itself, an array
        if fn_transform:
            return (
                numpy.array(cell_header                           , dtype=str),
                numpy.array(list(map(fn_transform, cell_metadata)), dtype=str)
            )

        return (
            numpy.array(cell_header  , dtype=str),
            numpy.array(cell_metadata, dtype=str)
        )


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
            return open(uncompressed_path, 'rb')

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

        def closure_cell_id_from_barcode(cell_record):
            """
            To get 'Cell ID' we prefix the cell barcode with '<sample_name>_' and we remove the
            '-1' from the end.
            """

            return '{}_{}'.format(sample_name, cell_record[0].rstrip('-1'))

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
        with cls.open_matrix_from_dir(path_to_mtx_dir) as matrix_handle:
            matrix_data = scipy.io.mmread(matrix_handle)

        gene_header, gene_metadata = GeneMetadataParser.from_handle(
            cls.open_genes_from_dir(path_to_mtx_dir),
            has_header=has_header
        )

        cell_header, cell_metadata = CellMetadataParser.from_handle(
            cls.open_cells_from_dir(path_to_mtx_dir),
            has_header=has_header
        )

        return (
            GeneExpression(
                expression_matrix=matrix_data.tocsc(),
                genes=gene_metadata[:,0],
                cells=cell_metadata[:,0]
            ),
            Annotation(
                headers=gene_header,
                annotations=gene_metadata
            ),
            Annotation(
                headers=cell_header,
                annotations=cell_metadata
            ),
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

        gene_header, gene_metadata = GeneMetadataParser.from_handle(
            cls.open_genes_from_dir(path_to_mtx_dir),
            has_header=has_header
        )

        barcode_header, barcode_metadata = CellMetadataParser.from_handle(
            cls.open_barcodes_from_dir(path_to_mtx_dir),
            has_header=has_header,
            fn_transform=closure_transform_barcode
        )

        return (
            GeneExpression(
                expression_matrix=matrix_data.tocsc(),
                genes=gene_metadata[:,0],
                cells=barcode_metadata[:,0]
            ),
            Annotation(
                headers=gene_header,
                annotations=gene_metadata
            ),
            Annotation(
                headers=barcode_header,
                annotations=barcode_metadata
            ),
        )
