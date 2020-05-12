#!/usr/bin/fish

set batch_sizes 100 1000 10000

#set base_search_dir /mnt/sdb/data/simulated-arrow-15/"batch-$batch_size"
set base_search_dir "/mnt/sdb/data/cost-model"
set pool_name       "cost-model"

for batch_size in $batch_sizes
    set search_dir "$base_search_dir/batch-$batch_size"

    for workload_dir in (find "$search_dir" -mindepth 1 -maxdepth 1 -type d)
        set num_files (ls -1 $workload_dir | wc -l)
        set table_name (basename $workload_dir)

        echo "$table_name-$batch_size ($num_files)"
        #echo "$ceph_opts"

        toolbox/rados-write --config-file $HOME/cluster/ceph.conf --input-dir $workload_dir --data-format arrow --use-skyhook-wrapper --pool $pool_name --table $table_name-$batch_size
    end
end
