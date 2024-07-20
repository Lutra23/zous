#!/bin/bash

# 帮助信息函数
print_usage() {
    echo "用法: $0"
}

# 检查输入参数
if [ "$#" -ne 0 ]; then
    echo "错误: 不应该有任何参数"
    print_usage
    exit 1
fi

# 寻找当前文件夹中的 gesp 文件
gesp_files=$(find . -maxdepth 1 -type f -name '*.gesp')

if [ -z "$gesp_files" ]; then
    echo "错误: 找不到 gesp 文件"
    exit 1
fi

# 显示可选的 gesp 文件列表
echo "找到以下 gesp 文件:"
count=0
for gesp_file in $gesp_files; do
    count=$((count+1))
    echo "$count. $gesp_file"
done

# 获取用户选择的文件序号
read -p "请选择要处理的文件 (输入序号): " choice

# 验证用户输入
if ! [[ "$choice" =~ ^[0-9]+$ ]]; then
    echo "错误: 无效的输入"
    exit 1
fi

if [ "$choice" -lt 1 ] || [ "$choice" -gt "$count" ]; then
    echo "错误: 选择的序号超出范围"
    exit 1
fi

# 获取用户选择的文件
selected_file=$(echo "$gesp_files" | sed -n "${choice}p")

# 提取文件名
mol2_name=$(basename "$selected_file" .gesp)

# 获取电荷信息
read -p "请输入原子的电荷信息: " charge

# 生成 mol2 文件
antechamber -i "$mol2_name.gesp" -fi gesp -o "$mol2_name.mol2" -fo mol2 -c resp -rn "$mol2_name" -eq 2

# 检查 mol2 文件是否生成成功
if [ ! -f "$mol2_name.mol2" ]; then
    echo "错误: 无法生成 mol2 文件: $mol2_name.mol2"
    exit 1
fi

# 生成 frcmod 文件
parmchk2 -i "$mol2_name.mol2" -f mol2 -o "$mol2_name.frcmod" -a "$charge"

# 检查 frcmod 文件是否生成成功
if [ ! -f "$mol2_name.frcmod" ]; then
    echo "错误: 无法生成 frcmod 文件: $mol2_name.frcmod"
    exit 1
fi

echo "已生成并移动文件到当前目录路径"
