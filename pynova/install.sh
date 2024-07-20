#!/bin/bash

# 函数：遍历目录中的所有文件并给予执行权限
traverse_and_chmod() {
    local dir="$1"
    # 遍历目录中的所有文件
    find "$dir" -type f -exec chmod +x {} \;
    echo "已给予目录 '$dir' 中的所有文件执行权限"
}

# 主程序
main() {

    # 遍历当前目录及子目录中的所有文件并给予执行权限
    traverse_and_chmod .

    # 执行其他必要的操作
    source /export/home/zoushuai/ben/name.sh
}

# 执行主程序
main

