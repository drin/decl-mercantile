#!/usr/bin/fish

set batch_sizes 100 1000 10000

for batch_size in $batch_sizes
    for workload_dir in (find /mnt/sdb/data/simulated-arrow-15/"batch-$batch_size" -mindepth 1 -maxdepth 1 -type d)
        set num_files (ls -1 $workload_dir | wc -l)
        set table_name (basename $workload_dir)

        # echo "$table_name ($num_files)"
        toolbox/rados-write --config-file $HOME/cluster/ceph.conf --input-dir $workload_dir  --data-format arrow --use-skyhook-wrapper --pool "simulated" --table "$table_name-$batch_size"
    end
end
