import numpy
import scipy.io

# ------------------------------
# file format parsers
class MatrixParser(object):

    __slots__ = ()

    @classmethod
    def from_path(cls, path_to_mtx):
        """
        Uses scipy's scipy.io.mmread() function to read the sparse matrix in MTX format.

        :path_to_mtx_root: Path to matrix data, in MTX format, with extensions: '.mtx' or
        '.mtx.gz'.
        """

        return scipy.io.mmread(path_to_mtx)


class GeneListParser(object):

    __slots__ = ()

    @classmethod
    def from_path(cls, path_to_genelist, has_header=True, delim='\t'):
        """
        This function assumes a two-column file where the first column is identical to the second. I
        don't know what the columns are, but I assume it is <gene symbol, gene alias>.
        """

        with open(path_to_genelist, 'r') as genelist_handle:
            # peel header, if necessary
            if has_header: next(genelist_handle)

            # hardcoded for census of immune cells for nwo
            # gid, sym, gene_type, seq_name, start_pos, end_pos, is_gene
            gene_list = [
                numpy.array(gene_line.strip().split(delim), dtype=str)
                for gene_line in genelist_handle
            ]

        return numpy.array(gene_list, dtype=str)

class CellIDParser(object):

    __slots__ = ()

    @classmethod
    def from_path(cls, path_to_cell_ids, has_header=True, delim='\t'):
        """
        This function assumes a one-column file containing only single-cell cell IDs. Currently, no
        validation of the file contents is done.
        """

        with open(path_to_cell_ids, 'r') as barcode_handle:
            # peel header, if necessary
            if has_header:
                metadata_columns = numpy.array(
                     next(barcode_handle).strip()
                                          .split(delim)
                    ,dtype=str
                )

            # eventually migrate to a dynamic access such as this
            # target_columns = numpy.where(metadata_columns in ['barcode',])

            # 'barcode' column has a 0-based index of 7
            barcode_col_ndx = 7
            cell_ids = [
                barcode_line.strip().split(delim)[barcode_col_ndx]
                for barcode_line in barcode_handle
            ]

        # return gene list as a numpy array of strings
        return numpy.array(cell_ids, dtype=str)

    @classmethod
    def from_moca(cls, path_to_barcodes, has_header=False, delim='\t'):
        pass


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
