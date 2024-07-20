#!/bin/bash

# 获取脚本所在目录的路径
folder_path="$(dirname "$BASH_SOURCE")"

# 遍历文件夹中的每个文件
for file in "$folder_path"/*; do
    # 检查是否为普通文件
    if [ -f "$file" ]; then
        # 创建别名
        alias_name=$(basename "$file")
        alias_name=${alias_name%.*}  # 去除文件扩展名
        alias "$alias_name"="$file"
    fi
done
