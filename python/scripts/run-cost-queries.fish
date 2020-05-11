set obj_count 10

# from objects with 100 columns
while test $obj_count -le 20
    echo "Querying $obj_count objects..."
    toolbox/query-benchmark --pool cost-model --num-objs $obj_count --query-str 'SELECT * FROM 1M-Immune-Cells-100'

    set obj_count (math $obj_count + 10)
end

# from objects with 1000 columns
while test $obj_count -le 20
    echo "Querying $obj_count objects..."
    toolbox/query-benchmark --pool cost-model --num-objs $obj_count --query-str 'SELECT * FROM 1M-Immune-Cells-1000'

    set obj_count (math $obj_count + 10)
end

#toolbox/query-benchmark --pool cost-model --num-objs 8697 --query-str 'SELECT * FROM 1M-Immune-Cells-100'
#toolbox/query-benchmark --pool cost-model --num-objs 866  --query-str 'SELECT * FROM 1M-Immune-Cells-1000'
#toolbox/query-benchmark --pool cost-model --num-objs 81  --query-str 'SELECT * FROM 1M-Immune-Cells-10000'
