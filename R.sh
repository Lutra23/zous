#!/bin/bash

# 设置颜色常量
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
BOLD_YELLOW='\033[1;33m' # 加粗的黄色
NC='\033[0m' # No Color

while true; do
    clear
	echo -e "${BLUE}######################################################################${NC}"
    echo -e "${BLUE}#      #####################################################         #${NC}"
    echo -e "${BLUE}#      #${NC}${BOLD_YELLOW} ██████╗    ███╗   ███╗    ███╗   ███╗  ███╗   ███╗${NC}${BLUE}#         #${NC}"
    echo -e "${BLUE}#      #${NC}${BOLD_YELLOW}██╔═══██╗   ████╗ ████║    ████╗ ████║  ████╗ ████║${NC}${BLUE}#         #${NC}"
    echo -e "${BLUE}#      #${NC}${BOLD_YELLOW}██║   ██║   ██╔████╔██║    ██╔████╔██║  ██╔████╔██║${NC}${BLUE}#         #${NC}"
    echo -e "${BLUE}#      #${NC}${BOLD_YELLOW}██║  ███║   ██║╚██╔╝██║    ██║╚██╔╝██║  ██║╚██╔╝██║${NC}${BLUE}#         #${NC}"
    echo -e "${BLUE}#      #${NC}${BOLD_YELLOW}╚████████║  ██║ ╚═╝ ██║    ██║ ╚═╝ ██║  ██║ ╚═╝ ██║${NC}${BLUE}#         #${NC}"
    echo -e "${BLUE}#      #${NC}${BOLD_YELLOW} ╚═════╝██║ ╚═╝     ╚═╝    ╚═╝     ╚═╝  ╚═╝     ╚═╝${NC}${BLUE}#         #${NC}"
    echo -e "${BLUE}######################################################################${NC}"
    echo -e "${RED}欢迎进入 QMMM 模拟操作系统!                                          #${NC}" 
    echo -e "${RED}                                                                     #${NC}"                                                                                     
    echo -e "${GREEN}请选择操作：                                                         #${NC}"
    echo -e "${BLUE}######################################################################${NC}"
    echo -e "${CYAN}1: 生成.com文件                                                      #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${CYAN}2: 生成mol2 frcmod文件                                               #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${CYAN}3: 去氢、加氢、合并、top crd合并为pdb                                #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${CYAN}4: 一键tleap                                                         #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${CYAN}5: 一键MD                                                            #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${CYAN}6: 取帧                                                              #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${CYAN}7: QM/MM 作业提交系统                                                #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${YELLOW}按Enter退出                                                          #${NC}"
    echo -e "${CYAN}                                                                     #${NC}"
    echo -e "${BLUE}######################################################################${NC}"
    read -r -p $'\033[1;33m请输入选项 (1-7): \033[0m' option


    case $option in
        1) sh $(dirname "$0")/runcom.sh ;;
        2) sh $(dirname "$0")/gresp.sh ;;
        3) sh $(dirname "$0")/pdb4amber.sh ;;
        4) sh $(dirname "$0")/leap.sh ;;
        5) sh $(dirname "$0")/MD.sh ;;
        6) sh $(dirname "$0")/frame.sh ;;
        7) python "$(dirname "$0")/main.py" ;;
        8 | "") 
            echo -e "${RED}退出脚本${NC}"
            break
            ;;
        *) echo -e "${RED}无效的选项，请重新输入${NC}" ;;
    esac

    read -n 1 -s -r -p "${GREEN}按任意键继续...${NC}"
    echo
done