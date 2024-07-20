import os

def xyz_to_gjf(xyz_filename, gjf_filename, title="Molecule", method="B3LYP", basis="6-31G(d)", charge=0, multiplicity=1):
    """
    将XYZ文件转换为Gaussian输入文件（.gjf）。
    
    :param xyz_filename: 输入XYZ文件名。
    :param gjf_filename: 输出GJF文件名。
    :param title: Gaussian输入文件的标题。
    :param method: 计算方法。
    :param basis: 基组。
    :param charge: 分子电荷。
    :param multiplicity: 自旋多重性。
    """
    with open(xyz_filename, 'r') as xyz_file:
        lines = xyz_file.readlines()

    # 跳过前两行，提取原子坐标行
    atom_lines = lines[2:]

    with open(gjf_filename, 'w') as gjf_file:
        # 写入路由部分
        gjf_file.write(f"%chk={gjf_filename[:-4]}.chk\n")
        gjf_file.write(f"# {method}/{basis} Opt\n\n")

        # 写入标题
        gjf_file.write(f"{title}\n\n")

        # 写入电荷和自旋多重性
        gjf_file.write(f"{charge} {multiplicity}\n")

        # 写入原子坐标
        for line in atom_lines:
            gjf_file.write(line)

        # 在文件末尾添加一个空行
        gjf_file.write("\n")

def batch_convert_xyz_to_gjf(folder_path, title="Molecule", method="B3LYP", basis="6-31G(d)", charge=0, multiplicity=1):
    """
    批量将文件夹中的所有XYZ文件转换为Gaussian输入文件（.gjf）。
    
    :param folder_path: 包含XYZ文件的文件夹路径。
    :param title: Gaussian输入文件的标题。
    :param method: 计算方法。
    :param basis: 基组。
    :param charge: 分子电荷。
    :param multiplicity: 自旋多重性。
    """
    for filename in os.listdir(folder_path):
        if filename.endswith('.xyz'):
            xyz_filepath = os.path.join(folder_path, filename)
            gjf_filename = filename.replace('.xyz', '.gjf')
            gjf_filepath = os.path.join(folder_path, gjf_filename)
            xyz_to_gjf(xyz_filepath, gjf_filepath, title, method, basis, charge, multiplicity)
            print(f'GJF文件已保存为 {gjf_filepath}')

def match_convert_xyz_to_gjf(folder_path):
    remove_files = [f for f in os.listdir(folder_path) if f.endswith('.xyz')]
    for file in remove_files:
        os.remove(os.path.join(folder_path, file))

    

if __name__ == '__main__':
    # 示例：指定包含XYZ文件的文件夹路径
    folder_path = r'C:\Users\szou5\Desktop\nn ns\XYZ'  # 替换为实际的文件夹路径
    batch_convert_xyz_to_gjf(folder_path)
    match_convert_xyz_to_gjf(folder_path)