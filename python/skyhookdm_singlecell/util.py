import argparse


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
            ,dest='input_metadata'
            ,type=str
            ,required=required
            ,help=(help_str or 'Path to file containing metadata for this program to process')
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

    def add_analysis_arg(self, required=False, help_str=''):
        self._arg_parser.add_argument(
             '--should-analyze'
            ,dest='should_analyze'
            ,action='store_true'
            ,required=required
            ,help=(help_str or 'Flag representing whether to do runtime analysis')
        )

        return self

    def parse_args(self):
        return self._arg_parser.parse_known_args()
