#!/bin/bash

# 定义物种数组
species=('niu' 'ji' 'guoying' 'dashu' 'ma' 'zhu' 'xianchong' 'monkey' 'banmayu' 'mianyang' 'quan' 'human' 'mouse')

# 遍历每个物种
for specie in "${species[@]}"; do
    echo "Processing specie: $specie"

    # 输入目录路径
    directory="/input/${specie}"

    # 创建输出目录（如果不存在的话）
    output_directory="/output/${specie}"
    mkdir -p "$output_directory"

    # 遍历目录下所有的CSV文件
    file_index=0
    find "$directory" -type f -name "*.csv" -print0 | while IFS= read -r -d '' csv_file; do
        ((file_index++))
        echo "Processing $file_index file: $csv_file"

        # 使用awk处理CSV文件，统计每一列零值的个数
        output=$(awk -F',' '{
            for (i=1; i<=NF; i++) {
                if ($i == 0) {
                    count[i]++
                }
            }
        } END {
            for (i=1; i<=NF; i++) {
                print (count[i] ? count[i] : 0)
            }
        }' "$csv_file")

        # 获取文件行数
        line_count=$(wc -l < "$csv_file")

        # 构建CSV文件内容
        csv_content=""
        i=1
        for zero_count in $output; do
            csv_content+="$i,$zero_count\n"
            ((i++))
        done

        # 构建输出文件名
        file_name=$(basename "$csv_file" .csv)
        output_file="${output_directory}/${file_name}_${line_count}_zero_count.csv"

        # 保存结果到CSV文件
        echo -e "$csv_content" > "$output_file"
        echo "Zero count for each column has been saved to $output_file."
    done
done
