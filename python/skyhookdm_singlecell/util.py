import sys
import logging
import argparse

import importlib


# ------------------------------
# Setup logger, primarily for optional imports
module_logger = logging.getLogger('util')
module_logger.addHandler(logging.StreamHandler(sys.stdout))
module_logger.setLevel(logging.DEBUG)


# ------------------------------
# convenience functions
def normalize_str(str_or_bytes, byte_encoding='utf-8'):
    if type(str_or_bytes) is bytes:
        return str_or_bytes.decode(byte_encoding)

    return str_or_bytes


def batched_indices(element_count, batch_size):
    for batch_id, batch_start in enumerate(range(0, element_count, batch_size), start=1):

        batch_end = batch_id * batch_size
        if batch_end > element_count:
            batch_end  = element_count

        yield batch_id, batch_start, batch_end


def try_import(module_name, is_required=False, module_package=None):
    imported_module = None

    try:
        imported_module = importlib.import_module(module_name, package=module_package)

    except ModuleNotFoundError as import_err:
        module_logger.warn(f'Error importing module {module_name}:\n{import_err}')

    finally:
        return imported_module


# ------------------------------
# utility classes
class ArgparseBuilder(object):

    @classmethod
    def with_description(cls, program_descr='Default program'):
        return cls(argparse.ArgumentParser(description=program_descr))

    def __init__(self, arg_parser, **kwargs):
        super(ArgparseBuilder, self).__init__(**kwargs)

        self._arg_parser = arg_parser

    def add_input_dir_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--input-dir'
            ,dest='input_dir'
            ,type=str
            ,required=required
            ,help=(help_str or
                   'Path to directory containing input files for this program to process')
        )

        return self

    def add_config_file_arg(self, required=False, help_str='', default=''):
        self._arg_parser.add_argument(
             '--config-file'
            ,dest='config_file'
            ,type=str
            ,required=required
            ,help=(help_str or 'Path to config file for this program to process')
        )

        return self

    def add_input_file_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--input-file'
            ,dest='input_file'
            ,type=str
            ,required=required
            ,help=(help_str or 'Path to file for this program to process')
        )

        return self

    def add_input_metadata_file_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--input-metadata'
            ,dest='input_metadata_file'
            ,type=str
            ,required=required
            ,help=(help_str or 'Path to file containing metadata for this program to process')
        )

        return self

    def add_metadata_target_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--metadata-target'
            ,dest='metadata_target'
            ,type=str
            ,required=required
            ,help=(help_str or 'Type of data described by provided metadata: <cells | genes>')
        )

        return self

    def add_delimiter_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--field-delimiter'
            ,dest='field_delimiter'
            ,type=str
            ,required=required
            ,default='\t'
            ,help=(help_str or 'String to use as delimiter between fields of each record (line)')
        )

        return self

    def add_data_format_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--data-format'
            ,dest='data_format'
            ,type=str
            ,required=required
            ,help=(help_str or 'Serialization format for data')
        )

        return self

    def add_output_file_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--output-file'
            ,dest='output_file'
            ,type=str
            ,required=required
            ,help=(help_str or 'Path to file for this program to produce')
        )

        return self

    def add_output_file_format_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--output-file-format'
            ,dest='output_file_format'
            ,type=str
            ,required=required
            ,help=(help_str or 'File format to use when serializing output data')
        )

        return self

    def add_output_dir_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--output-dir'
            ,dest='output_dir'
            ,type=str
            ,required=required
            ,help=(help_str or 'Path to directory for output files to be written')
        )

        return self

    def add_batch_args(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--batch-size'
            ,dest='batch_size'
            ,type=int
            ,required=required
            ,help=(help_str or 'Amount of entities to include in each batch (aka window size)')
        )

        self._arg_parser.add_argument(
             '--batch-offset'
            ,dest='batch_offset'
            ,type=int
            ,required=required
            ,help=(help_str or 'Amount of entities to shift the batch start (aka stride)')
        )

        return self

    def add_skyhook_query_bin_args(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--query-bin'
            ,dest='path_to_query_bin'
            ,type=str
            ,default='/mnt/sdb/skyhookdm-ceph/build/bin/run-query'
            ,required=required
            ,help=(help_str or "Path to skyhook's run-query binary")
        )

        return self

    def add_query_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--query-str'
            ,dest='query_str'
            ,type=str
            ,required=required
            ,help=(help_str or 'Structured query expression')
        )

        return self

    def add_skyhook_obj_slice_args(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--start-obj'
            ,dest='start_obj'
            ,type=str
            ,default='0'
            ,required=False
            ,help=(help_str or 'First partition to do object scan from')
        )

        self._arg_parser.add_argument(
             '--num-objs'
            ,dest='num_objs'
            ,type=str
            ,default='1000'
            ,required=required
            ,help=(help_str or 'Number of objects (including start) to scan')
        )

        return self

    def add_flatbuffer_flag_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--use-skyhook-wrapper'
            ,dest='flag_use_wrapper'
            ,action='store_true'
            ,required=required
            ,help=(help_str or "Whether to use skyhook's flatbuffer wrapper. <True | False>")
        )

        return self

    def add_has_header_flag_arg(self, required=False, help_str=''):
        default_help_msg = (
            "Specifies that cells and genes list files have header lines to skip. "
            "<True | False>"
        )

        self._arg_parser.add_argument(
             '--data-files-have-header'
            ,dest='flag_has_header'
            ,action='store_true'
            ,required=required
            ,help=(help_str or default_help_msg)
        )

        return self

    def add_prepend_header_flag_arg(self, required=False, help_str=''):
        default_help_msg = (
            "Specifies that a header for cells or genes annotation files should be inserted. "
            "<True | False>"
        )

        self._arg_parser.add_argument(
             '--prepend-metadata-header'
            ,dest='flag_add_header'
            ,action='store_true'
            ,required=required
            ,help=(help_str or default_help_msg)
        )

        return self

    def add_analysis_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--should-analyze'
            ,dest='should_analyze'
            ,action='store_true'
            ,required=required
            ,help=(help_str or 'Flag representing whether to do runtime analysis')
        )

        return self

    def add_ceph_pool_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--pool'
            ,dest='ceph_pool'
            ,type=str
            ,required=required
            ,help=(help_str or 'Name of Ceph pool to use.')
        )

        return self

    def add_skyhook_table_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--table'
            ,dest='skyhook_table'
            ,type=str
            ,required=required
            ,help=(help_str or 'Name of table to be used in SkyhookDM naming scheme')
        )

        return self

    def parse_args(self):
        return self._arg_parser.parse_known_args()
