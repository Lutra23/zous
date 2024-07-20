import os

def replace_title_and_remove_connectivity(input_file, new_title):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # 检查文件开头是否包含 %nprocshared= 和 %mem=58GB
    if not (lines[0].startswith("%nprocshared=") and lines[1].startswith("%mem=58GB")):
        # 如果缺少，则在文件开头添加这两行
        lines.insert(0, "%nprocshared=28\n")
        lines.insert(1, "%mem=58GB\n")

    # 遍历文件内容，替换标题
    for i, line in enumerate(lines):
        if line.strip().startswith("#"):
            lines[i] = f"# {new_title}\n"

    # 将处理后的内容写回文件
    with open(input_file, 'w') as f:
        f.writelines(lines)

def process_files_in_folder(folder_path, new_title):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gjf"):
            file_path = os.path.join(folder_path, file_name)
            replace_title_and_remove_connectivity(file_path, new_title)

folder_path = r"C:\Users\zous\Desktop\gjd"  # 替换为你的文件夹路径
new_title = "opt freq def2svp m062x scrf=(iefpcm,solvent=ccl4)"  # 替换为你想要的新标题
process_files_in_folder(folder_path, new_title)
