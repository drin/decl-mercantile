search_dir="/public/groups/xhca/data/single-cell/bpa/ebi"
for path_to_matrix_mtx in $(find $search_dir -name 'matrix.mtx' -type f); do
    input_dir=$(dirname "${path_to_matrix_mtx}")

    pushd "${input_dir}"

    gzip "matrix.mtx"
    gzip "cells.tsv"
    gzip "genes.tsv"

    popd
done
