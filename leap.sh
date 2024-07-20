#!/bin/bash

echo "请选择要运行的操作："
echo "1. 修改pdb复合物原子名称与mol2文件一致"
echo "2. 转换PDB文件为Amber格式"
echo "3. 生成tleap.in"
echo "4. 运行tleap.in"
echo "5. 二次tleap"
echo "6. 退出"

read -p "请输入选项编号：" choice

case $choice in
    1)
        echo "修改pdb复合物原子名称与mol2文件一致："
        python "$(dirname "$0")/amber-MD/pdb.py"
        ;;
    3)
        echo "生成tleap.in："
        python "$(dirname "$0")/amber-MD/leap.py"
        ;;

    4)
        echo "运行tleap.in："
        python "$(dirname "$0")/amber-MD/topcrd.py"
        ;;
    2)
        echo "查找PDB文件并转换为Amber格式："
        pdb_file=$(find . -maxdepth 1 -type f -name "*.pdb" | head -n 1)
        if [ -z "$pdb_file" ]; then
            echo "未找到PDB文件"
            exit 1
        fi
        pdb4amber -i "$pdb_file" -o "${pdb_file%.pdb}.amber.pdb"
        ;;
    5)
        echo "二次tleap："
        python "$(dirname "$0")/amber-MD/二次leap.py"
        ;;
    6)
        echo "退出脚本"
        exit 0
        ;;
    *)
        echo "无效选项，请重新输入"
        ;;
esac
