#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include "arrow/api.h"
#include "arrow/io/memory.h"
#include "arrow/ipc/reader.h"
#include "flatbuffers.h"

#include "fb_meta_generated.h"

namespace opt_parser = boost::program_options;

size_t size_in_bytes(std::ifstream *input_file_stream) {
    /*
     * Convenience function that takes an opened input file stream and:
     * - Seeks to the end
     * - Asks for the position (byte-based)
     * - Seeks to the beginning (for follow-up functions)
     */

    if (not input_file_stream->is_open()) { return 0; }

    input_file_stream->seekg(0, std::ios::end);
    size_t byte_count = input_file_stream->tellg();
    input_file_stream->seekg(0, std::ios::beg);

    return byte_count;
}

char* bytebuffer_from_filepath(std::string path_to_input) {
    /**
     * Opens file, <path_to_input>, and reads data into a char array (for use with arrow
     * functions). The opened file is closed when the input file stream (std::ifstream) goes out of
     * scope and is deconstructed.
     */

    std::cout << "[INFO] Parsing file: " << path_to_input << std::endl;

    std::ifstream input_file_stream(path_to_input, std::ios::in | std::ios::binary);
    if (not input_file_stream.is_open()) {
        std::cerr << "Unable to open file: " << path_to_input << std::endl;
        return NULL;
    }

    // Prep for reading file data into the char buffer
    size_t input_byte_count = size_in_bytes(&input_file_stream);
    char *data              = new char[input_byte_count];

    input_file_stream.read(data, input_byte_count);

    return data;
}

int table_from_record_reader(std::shared_ptr<arrow::ipc::RecordBatchReader> &record_reader,
                             std::shared_ptr<arrow::Table>                  &arrow_table) {

    std::vector<std::shared_ptr<arrow::RecordBatch>> accumulated_records;
    std::shared_ptr<arrow::RecordBatch>              current_record;

    // Iterate over RecordBatch objects using `record_reader` into an arrow Table
    while (true) {
        record_reader->ReadNext(&current_record);
        if (current_record == nullptr) { break; }

        accumulated_records.push_back(current_record);
    }

    arrow::Table::FromRecordBatches(accumulated_records, &arrow_table);

    return 0;
}

void parse_skyhook_format(const char *file_data) {
    std::shared_ptr<arrow::Buffer>                 arrow_buffer;
    std::shared_ptr<arrow::Table>                  arrow_table;
    std::shared_ptr<arrow::ipc::RecordBatchReader> arrow_record_reader;

    auto flatbuf_wrapper = Tables::GetFB_Meta(file_data);

    // Temporary Debug statements of some flatbuffer fields
    std::cout << "[FLATBUF] Blob Format: " << flatbuf_wrapper->blob_format()
              << std::endl
              << "[FLATBUF] Blob Size: "   << flatbuf_wrapper->blob_size()
              << std::endl;

    // Access the data from the flatbuffer metadata wrapper
    const flatbuffers::Vector<uint8_t> *stored_data_blob = flatbuf_wrapper->blob_data();

    // Create an arrow Buffer from a raw data string (constructed from raw storage blob)
    const std::string raw_buffer_data = std::string(
        reinterpret_cast<const char *>(stored_data_blob->data()),
        flatbuf_wrapper->blob_size()
    );
    arrow::Buffer::FromString(raw_buffer_data, &arrow_buffer);

    // Prepare a Buffer Reader (input stream) and RecordBatch Reader to stream over the arrow Buffer
    const std::shared_ptr<arrow::io::BufferReader> buffer_reader = std::make_shared<arrow::io::BufferReader>(arrow_buffer);
    arrow::ipc::RecordBatchStreamReader::Open(buffer_reader, &arrow_record_reader);

    // Validate that all is well by first grabbing the schema and metadata
    auto schema   = arrow_record_reader->schema();
    auto metadata = schema->metadata();

    std::cout << "[ARROW] Table Schema: "  << schema->ToString() << std::endl;
    // std::cout << "[ARROW Table Metadata: " << metadata << std::endl;

    // Read from `arrow_record_reader` and parse into `arrow_table`
    int parse_status = table_from_record_reader(arrow_record_reader, arrow_table);
    if (parse_status) {
        std::cerr << "[Error] Error when parsing arrow buffer into Table: "
                  << parse_status
                  << std::endl;
    }

    std::cout << "[ARROW] Table |> Row Count: " << arrow_table->num_rows()
              << std::endl;

    std::string column_name = "007e8251cae22f3198546f70ef52e91a";
    auto table_column = arrow_table->GetColumnByName(column_name);
    auto column_chunk = table_column->chunk(0);

    std::cout << "[ARROW] Table |> Column: "      << column_name
              << std::endl
              << "[ARROW] Table |> Chunk Count: " << table_column->num_chunks()
              << std::endl
              << "[ARROW] Chunk |> Data: "        << column_chunk->ToString()
              << std::endl;
}

int main (int argc, char **argv) {
    // Variables for parser to populate
    std::string path_to_input;

    opt_parser::variables_map       parsed_args;
    opt_parser::options_description parser_opts("Parser Options");

    // Define CLI options
    parser_opts.add_options()
        (
            "help,h",
            "show help message"
        )
        (
            "input-file",
            opt_parser::value<std::string>(&path_to_input)->default_value(
                "/raiddata/data/skyhook/example.skyhook-blob"
            ),
            "Path to skyhook file to parse. "
            "(default: '/raiddata/data/skyhook/example.skyhook-blob')"
        )
    ;

    // Parse CLI args into variables_map
    opt_parser::store(opt_parser::parse_command_line(argc, argv, parser_opts), parsed_args);

    // If help is specified as an option
    if (parsed_args.count("help")) {
        // print the defined options
        std::cout << parser_opts << std::endl;
        return 1;
    }

    opt_parser::notify(parsed_args);

    const char *file_data = bytebuffer_from_filepath(path_to_input);

    if (file_data == NULL) {
        std::cerr << "File '" << path_to_input << "' could not be parsed"
                  << std::endl;
    }

    // Pass file data to be parsed as a skyhook flatbuffer metadata wrapper
    parse_skyhook_format(file_data);

    return 0;
}
