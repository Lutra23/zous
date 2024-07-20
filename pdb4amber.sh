#!/bin/bash

# 设置颜色常量
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

while true; do
    clear
    echo -e "${GREEN}请选择要执行的命令：${NC}"
    echo -e "${YELLOW}1: 使用reduce命令处理 .pdb 文件${NC}"
    echo -e "${YELLOW}2: 使用cat命令合并 .pdb 文件${NC}"
    echo -e "${YELLOW}3: 使用ambpdb命令转换 .pdb 文件${NC}"
    echo -e "${YELLOW}4: 执行去氢命令${NC}"
    echo -e "${YELLOW}5: 退出${NC}"

    read -p "请输入选项 (1-5): " option

    case $option in
        1) 
            echo -e "${GREEN}请选择要处理的 .pdb 文件:${NC}"
            # 查找当前目录中的所有 .pdb 文件
            pdb_files=$(find . -maxdepth 1 -type f -name "*.pdb")

            # 检查是否存在 .pdb 文件
            if [ -z "$pdb_files" ]; then
                echo -e "${RED}当前文件夹中找不到任何 .pdb 文件${NC}"
                continue
            fi

            # 创建一个数组来保存 .pdb 文件列表
            pdb_list=()
            index=1

            # 将 .pdb 文件添加到列表中
            for file in $pdb_files; do
                echo "$index: $file"
                pdb_list+=("$file")
                ((index++))
            done

            # 让用户选择要处理的 .pdb 文件
            read -p "请输入选项 (1-${#pdb_list[@]}): " pdb_option

            # 检查用户选择的选项是否有效
            if [ "$pdb_option" -ge 1 ] && [ "$pdb_option" -le "${#pdb_list[@]}" ]; then
                selected_pdb="${pdb_list[$((pdb_option-1))]}"
                echo -e "${BLUE}正在执行reduce命令...${NC}"
                reduce "$selected_pdb" > "${selected_pdb%.pdb}_H.pdb"
                echo -e "${GREEN}处理完成，生成文件: ${selected_pdb%.pdb}_H.pdb${NC}"
            else
                echo -e "${RED}无效的选项，请重新选择${NC}"
            fi
            ;;
        2)
            echo -e "${GREEN}请选择要合并的 .pdb 文件:${NC}"
            # 查找当前目录中的所有 .pdb 文件
            pdb_files=$(find . -maxdepth 1 -type f -name "*.pdb")

            # 检查是否存在 .pdb 文件
            if [ -z "$pdb_files" ]; then
                echo -e "${RED}当前文件夹中找不到任何 .pdb 文件${NC}"
                continue
            fi

            # 创建一个数组来保存选择的 .pdb 文件
            selected_files=()

            # 循环提示用户选择 .pdb 文件，直到输入空值为止
            while true; do
                # 打印 .pdb 文件列表
                echo -e "${YELLOW}当前已选择的 .pdb 文件:${NC}"
                for file in "${selected_files[@]}"; do
                    echo "  - $file"
                done
                echo

                echo -e "${YELLOW}请选择要合并的 .pdb 文件 (输入空值以结束选择):${NC}"
                select pdb_file in $pdb_files "取消"; do
                    if [ "$pdb_file" = "取消" ]; then
                        break
                    elif [ ! -z "$pdb_file" ]; then
                        selected_files+=("$pdb_file")
                        break
                    else
                        echo -e "${RED}无效的选项，请重新选择${NC}"
                    fi
                done
            done

            # 执行cat命令
            if [ ${#selected_files[@]} -gt 0 ]; then
                echo -e "${BLUE}正在执行cat命令...${NC}"
                cat "${selected_files[@]}" > concatenated.pdb
                echo -e "${GREEN}合并完成，生成文件: concatenated.pdb${NC}"
            else
                echo -e "${RED}未选择任何 .pdb 文件，取消操作${NC}"
            fi
            ;;
        3) 
            echo -e "${GREEN}请选择 .top 文件:${NC}"
            # 查找当前目录中的所有 .top 文件
            top_files=$(find . -maxdepth 1 -type f -name "*.prmtop")

            # 检查是否存在 .top 文件
            if [ -z "$top_files" ]; then
                echo -e "${RED}当前文件夹中找不到任何 .top 文件${NC}"
                continue
            fi

            # 打印 .top 文件列表
            select top_file in $top_files; do
                if [ ! -z "$top_file" ]; then
                    echo -e "${GREEN}已选择 .top 文件: ${top_file}${NC}"
                    break
                else
                    echo -e "${RED}无效的选项，请重新选择${NC}"
                fi
            done

            echo -e "${GREEN}请选择 .crd 文件:${NC}"
            # 查找当前目录中的所有 .crd 文件
            crd_files=$(find . -maxdepth 1 -type f -name "*.inpcrd")

            # 检查是否存在 .crd 文件
            if [ -z "$crd_files" ]; then
                echo -e "${RED}当前文件夹中找不到任何 .crd 文件${NC}"
                continue
            fi

            # 打印 .crd 文件列表
            select crd_file in $crd_files; do
                if [ ! -z "$crd_file" ]; then
                    echo -e "${GREEN}已选择 .crd 文件: ${crd_file}${NC}"
                    break
                else
                    echo -e "${RED}无效的选项，请重新选择${NC}"
                fi
            done

            # 执行ambpdb命令
            echo -e "${BLUE}正在执行ambpdb命令...${NC}"
            ambpdb -p "$top_file" -c "$crd_file" > converted.pdb
            echo -e "${GREEN}转换完成，生成文件: converted.pdb${NC}"
            ;;
        4) 
            echo -e "${GREEN}正在执行去氢命令...${NC}"
            # 查找当前目录中的所有 .pdb 文件
            pdb_files=$(find . -maxdepth 1 -type f -name "*.pdb")

            # 检查是否存在 .pdb 文件
            if [ -z "$pdb_files" ]; then
                echo -e "${RED}当前文件夹中找不到任何 .pdb 文件${NC}"
                continue
            fi

            # 创建一个数组来保存选择的 .pdb 文件
            pdb_list=()
            index=1

            # 将 .pdb 文件添加到列表中
            for file in $pdb_files; do
                echo "$index: $file"
                pdb_list+=("$file")
                ((index++))
            done

            # 让用户选择要处理的 .pdb 文件
            read -p "请输入选项 (1-${#pdb_list[@]}): " pdb_option

            # 检查用户选择的选项是否有效
            if [ "$pdb_option" -ge 1 ] && [ "$pdb_option" -le "${#pdb_list[@]}" ]; then
                selected_pdb="${pdb_list[$((pdb_option-1))]}"
                echo -e "${BLUE}正在执行reduce命令...${NC}"
                reduce "$selected_pdb" > "${selected_pdb%.pdb}_H.pdb"
                echo -e "${GREEN}去氢完成，生成文件: ${selected_pdb%.pdb}_H.pdb${NC}"
            else
                echo -e "${RED}无效的选项，请重新选择${NC}"
            fi
            ;;
        5 | "") 
            echo -e "${RED}退出脚本${NC}"
            break
            ;;
        *) 
            echo -e "${RED}无效的选项，请重新输入${NC}"
            ;;
    esac

    read -n 1 -s -r -p "${GREEN}按任意键继续...${NC}"
done
