import subprocess
import os
import shutil

def copy_files(source_folder, destination_folder):
    
    # 遍历源文件夹中的所有文件，并复制到目标文件夹
    for filename in os.listdir(source_folder):
        source_file_path = os.path.join(source_folder, filename)
        destination_file_path = os.path.join(destination_folder, filename)
        shutil.copy2(source_file_path, destination_file_path)

# 调用函数，将源文件夹中的所有文件复制到目标文件夹
source_folder = os.path.dirname(__file__) + "/scripts/函数/初始化/"
destination_folder = os.path.dirname(__file__) + "/scripts/函数/"
copy_files(source_folder, destination_folder)


def remove_tcl_files(folder_path):

    # 遍历目标文件夹中的所有文件，并删除所有tcl文件
    for filename in os.listdir(folder_path):
        if filename.endswith(".tcl"):
            file_path = os.path.join(folder_path, filename)
            os.remove(file_path)

# 调用函数，删除目标文件夹中的所有tcl文件
target_folder = os.path.dirname(__file__) + "/scripts/函数/"
remove_tcl_files(target_folder)


def run_script(script_name):
    script_path = os.path.dirname(__file__) + "/scripts/" + script_name
    subprocess.run(["python", script_path])

def main():
    print("========== 欢迎使用 QMMM 作业提交系统 ==========")
    print("请选择要运行的脚本：")
    print("1. QM/MM 小功能")
    print("2. QM/MM 泛函基组修改")
    print("3. QM/MM 分子优化")
    print("4. QM/MM 限制性优化")
    print("5. QM/MM scan 扫描")
    print("6. QM/MM 过渡态寻找 (prfo法)")
    print("7. QM/MM 过渡态寻找 (dimer法)")
    print("8. QM/MM 二次作业")
    print("9. 按 Enter 退出")

    choice = input("请输入数字选择 (1-9): ")
    if choice == "1":
        run_script("sub.py")
    elif choice == "2":
        run_script("basic.py")
    elif choice == "3":
        run_script("optrun.py")
    elif choice == "4":
        run_script("optfre.py")
    elif choice == "5":
        run_script("scanrun.py")
    elif choice == "6":
        run_script("tsrun.py")
    elif choice == "7":
        run_script("ts dimer.py")
    elif choice == "8":
        run_script("twice opt.py")
    elif choice == "9":
        print("感谢使用 QMMM 作业提交系统，再见！")
    else:
        print("无效的选择，请输入数字 1-9。")

if __name__ == "__main__":
    main()
