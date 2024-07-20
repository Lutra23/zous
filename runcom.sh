#!/bin/bash

# 执行第一个脚本
"$(dirname "$0")/amber-MD/com.sh"

# 检查第一个脚本是否成功执行
if [ $? -eq 0 ]; then
    echo "第一个脚本执行成功。"
    read -p "是否要继续执行第二个脚本？(Y/N): " choice
    case "$choice" in
        y|Y) 
            echo "正在执行第二个脚本..."
            python "$(dirname "$0")/amber-MD/convert_com.py"
            ;;
        n|N)
            echo "用户选择不继续执行第二个脚本。"
            ;;
        *)
            echo "无效的输入。"
            ;;
    esac
else
    echo "第一个脚本执行失败，无法执行第二个脚本。"
fi
