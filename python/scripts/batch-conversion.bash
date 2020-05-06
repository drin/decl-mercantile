#output_root_dir="/public/groups/xhca/data/outputs/arrow-15/batch-100.uint16"
#output_root_dir="/public/groups/xhca/data/outputs/simulated/batch-10000"
output_root_dir="/public/groups/xhca/data/outputs/cost-model"

#batch_sizes=(100 1000 10000)
batch_sizes=(100)

# ------------------------------
# For HCA data

#search_dir="/public/groups/xhca/data/single-cell/hca/human/Fetal-Maternal-Interface"
search_dir="/public/groups/xhca/data/single-cell/hca/human"
for batch_size in ${batch_sizes[@]}; do
    for path_to_matrix_mtx in $(find $search_dir -name 'matrix.mtx*' -type f); do
        input_dir=$(dirname "${path_to_matrix_mtx}")

        dataset_name=$(basename $(dirname "${input_dir}"))
        echo "${dataset_name}"

        output_dir="${output_root_dir}/batch-${batch_size}/${dataset_name}"

        mkdir -p "${output_dir}"

        conversion_opts="--batch-size ${batch_size} --output-file-format arrow --data-format arrow --data-files-have-header"
        toolbox/convert-file-format --input-dir "${input_dir}" --output-dir "${output_dir}" ${conversion_opts} 
    done
done


# ------------------------------
# For Ding/Bianca's data

##search_dir="/public/groups/xhca/data/single-cell/bpa/ebi"
#search_dir="/public/groups/xhca/workloads/simulated-data/datasets"
#for batch_size in ${batch_sizes[@]}; do
#    for path_to_matrix_mtx in $(find $search_dir -name 'matrix.mtx*' -type f); do
#        input_dir=$(dirname "${path_to_matrix_mtx}")
#
#        dataset_name=$(basename "${input_dir}")
#        echo "${dataset_name}"
#
#        output_dir="${output_root_dir}/batch-${batch_size}/${dataset_name}"
#
#        mkdir -p "${output_dir}"
#
#        conversion_opts="--batch-size ${batch_size} --output-file-format arrow --data-format arrow"
#        toolbox/convert-file-format --input-dir "${input_dir}" --output-dir "${output_dir}" ${conversion_opts} 
#        break
#    done
#done
