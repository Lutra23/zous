import os
import logging
from concurrent.futures import ProcessPoolExecutor

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
    try:
        with open(xyz_filename, 'r') as xyz_file:
            lines = xyz_file.readlines()

        atom_lines = lines[2:]

        with open(gjf_filename, 'w') as gjf_file:
            gjf_file.write(f"%chk={gjf_filename[:-4]}.chk\n")
            gjf_file.write(f"# {method}/{basis} Opt\n\n")
            gjf_file.write(f"{title}\n\n")
            gjf_file.write(f"{charge} {multiplicity}\n")

            for line in atom_lines:
                gjf_file.write(line)

            gjf_file.write("\n")
        logging.info(f'GJF文件已保存为 {gjf_filename}')
    except IOError as e:
        logging.error(f"文件操作失败: {e}")
    except Exception as e:
        logging.error(f"转换失败: {e}")

def create_new_folder(base_folder, folder_index):
    new_folder = os.path.join(base_folder, f'batch_{folder_index}')
    try:
        os.makedirs(new_folder, exist_ok=True)
        logging.info(f"创建文件夹: {new_folder}")
    except FileExistsError:
        logging.warning(f"文件夹 {new_folder} 已存在")
    except PermissionError:
        logging.error(f"没有权限创建文件夹 {new_folder}")
    except Exception as e:
        logging.error(f"创建文件夹失败: {e}")
    return new_folder

def process_file(args):
    xyz_filepath, output_folder, title, method, basis, charge, multiplicity = args
    filename = os.path.basename(xyz_filepath)
    gjf_filename = filename.replace('.xyz', '.gjf')
    gjf_filepath = os.path.join(output_folder, gjf_filename)
    xyz_to_gjf(xyz_filepath, gjf_filepath, title, method, basis, charge, multiplicity)

def batch_convert_xyz_to_gjf(folder_path, base_output_path, title="Molecule", method="B3LYP", basis="def2svp", charge=0, multiplicity=1, batch_size=10000):
    """
    批量将文件夹中的所有XYZ文件转换为Gaussian输入文件（.gjf），并每50000个文件创建一个新的文件夹。
    
    :param folder_path: 包含XYZ文件的文件夹路径。
    :param base_output_path: 基础输出文件夹路径。
    :param title: Gaussian输入文件的标题。
    :param method: 计算方法。
    :param basis: 基组。
    :param charge: 分子电荷。
    :param multiplicity: 自旋多重性。
    :param batch_size: 每个输出文件夹包含的文件数量。
    """
    try:
        xyz_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.xyz')]
        folder_index = 1
        output_folder = create_new_folder(base_output_path, folder_index)
        file_count = 0
        tasks = []

        with ProcessPoolExecutor() as executor:
            for xyz_filepath in xyz_files:
                if file_count >= batch_size:
                    folder_index += 1
                    output_folder = create_new_folder(base_output_path, folder_index)
                    file_count = 0

                task_args = (xyz_filepath, output_folder, title, method, basis, charge, multiplicity)
                tasks.append(task_args)
                file_count += 1

            executor.map(process_file, tasks)
    except Exception as e:
        logging.error(f"批处理转换失败: {e}")

if __name__ == '__main__':
    folder_path = r'C:\Users\szou5\Desktop\work\处理分子\xyz_files'  # 替换为实际的文件夹路径
    base_output_path = r'C:\Users\szou5\Desktop\work\处理分子\gjf'  # 替换为实际的基础输出文件夹路径
    batch_convert_xyz_to_gjf(folder_path, base_output_path)
