# core libraries
import sys

# dependencies
import pyarrow


# ------------------------------
# Module-level Variables
debug = True


# ------------------------------
# Functions
def write_data_as_arrow(path_to_outfile, batch_size, batch_offset, data_schema, data_recordbatch):
    batch_ndx = 0

    print('batch size: {}'.format(batch_size))
    print('batch offset: {}'.format(batch_offset))

    with open('{}'.format(path_to_outfile), 'wb') as output_handle:
        batch_writer = pyarrow.RecordBatchStreamWriter(output_handle, schema=data_schema)

        while batch_ndx + batch_offset < data_recordbatch.num_rows:

            if debug:
                sys.stdout.write('\rwriting batch: {:0d}'.format(batch_ndx))
                sys.stdout.flush()

            batch_writer.write_batch(data_recordbatch.slice(offset=batch_ndx, length=batch_size))
            batch_ndx += batch_offset

        # have to serialize the remainder batch (num_rows % batch_offset)
        if data_recordbatch.num_rows > batch_ndx:
            sys.stdout.write('\rwriting batch: {:0d}'.format(batch_ndx))
            sys.stdout.flush()

            batch_writer.write_batch(data_recordbatch.slice(offset=batch_ndx, length=batch_size))

        print('--- batches written')
