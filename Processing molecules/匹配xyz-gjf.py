import os
import shutil

def get_prefixes_from_xyz_files(xyz_folder_path):
    prefixes = set()
    for filename in os.listdir(xyz_folder_path):
        if filename.endswith('.xyz'):
            prefix = os.path.splitext(filename)[0]
            prefixes.add(prefix)
    return prefixes

def move_unmatched_gjf_files(gjf_folder, target_folder, xyz_prefixes):
    for filename in os.listdir(gjf_folder):
        if filename.endswith('.gjf'):
            prefix = os.path.splitext(filename)[0]
            if prefix not in xyz_prefixes:
                source_file_path = os.path.join(gjf_folder, filename)
                target_file_path = os.path.join(target_folder, filename)
                shutil.move(source_file_path, target_file_path)
                print(f"Moved '{filename}' to '{target_folder}'")

# 指定文件夹路径
source_folder = r"C:\Users\szou5\Desktop\work\gjf\reasonable"  # 替换为您的.gjf文件所在文件夹路径
xyz_folder = r"C:\Users\szou5\Desktop\work\处理分子\xyz"  # .xyz文件的文件夹路径
target_folder = r"C:\Users\szou5\Desktop\work\gjf\unmatched"  # 替换为目标文件夹路径

# 获取.xyz文件的前缀
xyz_prefixes = get_prefixes_from_xyz_files(source_folder)

# 移动没有匹配到的.gjf文件
move_unmatched_gjf_files(source_folder, target_folder, xyz_prefixes)

