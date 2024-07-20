import os
import re

def replace_title_and_handle_charge_spin(input_file, new_title, charge, spin):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 提取文件名（不包含扩展名）
    file_name_without_ext = os.path.splitext(os.path.basename(input_file))[0]
    chk_line = f"%chk={file_name_without_ext}.chk\n"

    # 检查文件开头是否包含 %nprocshared= 和 %mem=58GB
    if not (lines[0].startswith("%nprocshared=") and lines[1].startswith("%mem=58GB")):
        # 如果缺少，则在文件开头添加这两行
        lines.insert(0, "%nprocshared=28\n")
        lines.insert(1, "%mem=58GB\n")

    # 遍历文件内容，替换第一个以 # 开头的行作为新标题，并添加/替换 %chk= 行
    title_replaced = False
    chk_replaced = False
    charge_spin_line = f"{charge} {spin}\n"
    charge_spin_replaced = False
    molecule_section = False

    for i, line in enumerate(lines):
        if line.strip().startswith("#") and not title_replaced:
            lines[i] = f"# {new_title}\n"
            title_replaced = True
        elif line.strip().startswith("%chk="):
            lines[i] = chk_line
            chk_replaced = True
        elif molecule_section and re.match(r'^\s*-?\d+\s+-?\d+\s*$', line.strip()):
            lines[i] = charge_spin_line
            charge_spin_replaced = True
        elif line.strip() == "Molecule":
            molecule_section = True

    if not chk_replaced:
        # 如果没有找到 %chk= 行，则在文件开头添加
        lines.insert(0, chk_line)

    if not charge_spin_replaced:
        # 如果没有找到电荷和自旋行，则在适当位置添加
        for i, line in enumerate(lines):
            if line.strip() == "Molecule":
                lines.insert(i + 2, charge_spin_line)
                break

    # 将处理后的内容写回文件
    with open(input_file, 'w', encoding='utf-8') as f:
        f.writelines(lines)

def process_files_in_folder(folder_path, new_title, charge, spin):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gjf"):
            file_path = os.path.join(folder_path, file_name)
            replace_title_and_handle_charge_spin(file_path, new_title, charge, spin)

folder_path = r"C:\Users\szou5\Desktop\work\处理分子\gjf\batch\unreasonable"  # 替换为你的文件夹路径
new_title = "opt freq def2svp m062x scrf=(iefpcm,solvent=THF)"  # 替换为你想要的新标题
charge = 0  # 替换为你想要的电荷
spin = 1  # 替换为你想要的自旋
process_files_in_folder(folder_path, new_title, charge, spin)
