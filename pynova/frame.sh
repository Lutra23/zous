#!/bin/bash

while true; do
    # 交互式获取操作选项
    echo "请选择操作："
    echo "0: 提取特定帧的结构信息"
    echo "1: 提取特定区间的所有帧"
    echo "2: 将输入的 .mdcrd 文件转换成 .dcd 格式的输出文件"
    echo "3: 退出"
    read -p "请输入选项 (0, 1, 2 或 3): " option

    if [ $option -eq 0 ]; then
        # 提取特定帧的结构信息
        read -p "请输入要提取结构信息的帧数: " frame
        
  # 自动搜索当前文件夹内的拓扑文件
        prmtop_file=$(find "$PWD" -maxdepth 1 -type f -name '*.prmtop' -print -quit)

        if [ -z "$prmtop_file" ]; then
            echo "找不到 .prmtop 文件"
            continue
        fi

        # 获取当前文件夹中的所有 .dcd 文件
        dcd_files=($(find "$PWD" -maxdepth 1 -type f -name '*.dcd'))

        if [ ${#dcd_files[@]} -eq 0 ]; then
            echo "当前文件夹中找不到任何 .dcd 文件"
            continue
        fi

        # 显示找到的 .dcd 文件列表，供用户选择
        echo "找到的 .dcd 文件如下："
        for ((i=0; i<${#dcd_files[@]}; i++)); do
            echo "$i: ${dcd_files[i]}"
        done

        # 让用户选择要处理的 .dcd 文件
        read -p "请选择要处理的 .dcd 文件编号: " file_index

        if [ $file_index -ge 0 ] && [ $file_index -lt ${#dcd_files[@]} ]; then
            selected_file=${dcd_files[file_index]}
            # 提取指定帧的结构信息
            cpptraj $prmtop_file << EOF
            trajin wat_md_1.mdcrd $frame $frame
            trajout ${frame}.rst restartd
EOF
            ambpdb -p $prmtop_file -c ${frame}.rst > ${frame}.pdb

        else
            echo "无效的文件编号"
        fi
    elif [ $option -eq 1 ]; then
        # 提取特定区间的所有帧
        # 自动搜索当前文件夹内的拓扑文件
        prmtop_file=$(find "$PWD" -maxdepth 1 -type f -name '*.prmtop' -print -quit)

        if [ -z "$prmtop_file" ]; then
            echo "找不到 .prmtop 文件"
            continue
        fi

        # 获取当前文件夹中的所有 .dcd 文件
        dcd_files=($(find "$PWD" -maxdepth 1 -type f -name '*.dcd'))

        if [ ${#dcd_files[@]} -eq 0 ]; then
            echo "当前文件夹中找不到任何 .dcd 文件"
            continue
        fi

        # 显示找到的 .dcd 文件列表，供用户选择
        echo "找到的 .dcd 文件如下："
        for ((i=0; i<${#dcd_files[@]}; i++)); do
            echo "$i: ${dcd_files[i]}"
        done

        # 让用户选择要处理的 .dcd 文件
        read -p "请选择要处理的 .dcd 文件编号: " file_index

        if [ $file_index -ge 0 ] && [ $file_index -lt ${#dcd_files[@]} ]; then
            selected_file=${dcd_files[file_index]}

            # 交互式获取起始帧和结束帧
            read -p "请输入要提取的起始帧数: " start_frame
            read -p "请输入要提取的结束帧数: " end_frame

            # 提取指定区间内的所有帧
            cpptraj -p $prmtop_file << EOF
            trajin $selected_file $start_frame $end_frame
            trajout extracted_frames.dcd
EOF
            echo "提取特定区间 $start_frame 到 $end_frame 的所有帧完成，并保存为 extracted_frames.dcd"
        else
            echo "无效的文件编号"
        fi
    elif [ $option -eq 2 ]; then
        # 将输入的 .mdcrd 文件转换成 .dcd 格式的输出文件
        # 自动搜索当前文件夹中的拓扑文件
        prmtop_file=$(find "$PWD" -maxdepth 1 -type f -name '*.prmtop' -print -quit)

        if [ -z "$prmtop_file" ]; then
            echo "找不到 .prmtop 文件"
            continue
        fi

        # 获取当前文件夹中的所有 .mdcrd 文件
        mdcrd_files=($(find "$PWD" -maxdepth 1 -type f -name '*.mdcrd'))

        if [ ${#mdcrd_files[@]} -eq 0 ]; then
            echo "当前文件夹中找不到任何 .mdcrd 文件"
            continue
        fi

        # 显示找到的 .mdcrd 文件列表，供用户选择
        echo "找到的 .mdcrd 文件如下："
        for ((i=0; i<${#mdcrd_files[@]}; i++)); do
            echo "$i: ${mdcrd_files[i]}"
        done

        # 让用户选择要处理的 .mdcrd 文件
        read -p "请选择要转换的 .mdcrd 文件编号: " file_index

        if [ $file_index -ge 0 ] && [ $file_index -lt ${#mdcrd_files[@]} ]; then
            selected_file=${mdcrd_files[file_index]}

            # 自动生成输出文件名
            output_file="${selected_file%.*}.dcd"

            # 执行文件格式转换
            cpptraj -p $prmtop_file << EOF
            trajin $selected_file
            trajout $output_file
EOF
            echo "文件格式转换完成，已保存为 $output_file"
        else
            echo "无效的文件编号"
        fi
    elif [ $option -eq 3 ]; then
        # 退出
        echo "退出脚本"
        break
    else
        echo "无效的选项"
    fi
done
