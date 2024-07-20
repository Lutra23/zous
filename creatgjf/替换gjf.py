import os
import re

def replace_title_and_handle_charge_spin(input_file, new_title, charge, spin):
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        with open(input_file, 'r', encoding='GBK') as f:
            lines = f.readlines()

    file_name_without_ext = os.path.splitext(os.path.basename(input_file))[0]
    chk_line = f"%chk={file_name_without_ext}.chk\n"

    if not (lines[0].startswith("%nprocshared=") and lines[1].startswith("%mem=58GB")):
        lines.insert(0, "%nprocshared=28\n")
        lines.insert(1, "%mem=58GB\n")

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
        lines.insert(0, chk_line)

    if not charge_spin_replaced:
        for i, line in enumerate(lines):
            if line.strip() == "Molecule":
                lines.insert(i + 1, charge_spin_line)
                break

    with open(input_file, 'w', encoding='utf-8') as f:
        f.writelines(lines)

def process_files_in_folder(folder_path, new_title, charge, spin):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gjf"):
            file_path = os.path.join(folder_path, file_name)
            replace_title_and_handle_charge_spin(file_path, new_title, charge, spin)

# 以下路径、标题、电荷和自旋需要根据实际情况替换
folder_path = r"C:\Users\szou5\Desktop\nn ns\XYZ\unreasonable"  # 替换为你的文件夹路径
new_title = "opt freq def2svp m062x scrf=(iefpcm,solvent=THF)"  # 替换为你想要的新标题
charge = -1  # 替换为你想要的电荷
spin = 1  # 替换为你想要的自旋

process_files_in_folder(folder_path, new_title, charge, spin)
