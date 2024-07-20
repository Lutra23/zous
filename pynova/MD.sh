#!/bin/bash

# 脚本路径
update_simulation_steps_script=""$(dirname "$0")/amber-MD/run2.sh""
update_configuration_files_script=""$(dirname "$0")/amber-MD/run3.sh""
current_directory="$(dirname "$0")"
run_molecular_dynamics_script=""$(dirname "$0")/amber-MD/run.sh""

while true; do
    # 显示任务选择菜单
    echo "请选择要执行的任务:"
    echo "1. 更新模拟步数"
    echo "2. 更新残基以及约束力"
    echo "3. 运行分子动力学模拟"
    echo "4. 退出"

    # 读取用户选择
    read -p "请输入任务编号: " choice

    # 根据用户选择调用相应的脚本
    case $choice in
        1)
            bash "$update_simulation_steps_script"
            ;;
        2)
            bash "$update_configuration_files_script"
            ;;
        3)
            # 检查当前GPU是否有任务
            if nvidia-smi | grep "No running processes found" > /dev/null; then
                # 如果没有任务，提示用户继续
                read -p "当前GPU无任务，按 Enter 键继续运行，或按任意键后按Enter键取消..." continue_key
                if [ -z "$continue_key" ]; then
                    # 用户按下 Enter 键，继续运行
                    nohup bash "$run_molecular_dynamics_script" &
                else
                    # 用户按下了其他键，取消操作
                    echo "取消运行."
                fi
            else
                # 如果有任务，显示任务信息，并提示用户按 Enter 键上交任务
                echo "当前GPU有以下任务运行:"
                nvidia-smi
                read -p "按 Enter 键上交任务，或按任意键后按Enter键取消..." continue_key
                if [ -z "$continue_key" ]; then
                    # 用户按下 Enter 键，上交任务
                    echo "正在上交任务..."
                    # 这里执行上交任务的命令，此处仅为示例
                    echo "任务已上交."
                else
                    # 用户按下了其他键，取消操作
                    echo "取消上交任务."
                fi
            fi
            ;;
        4)
            echo "退出脚本."
            exit 0
            ;;
        *)
            echo "无效的选择，请重新输入."
            ;;
    esac
done
