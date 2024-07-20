import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from typing import List

# 定义一些常量以提高代码的可读性
DEFAULT_CHUNK_SIZE = 10000

def validate_path(path: str) -> str:
    """
    验证路径的有效性并返回清理后的路径。
    这是一个简单的示例，实际使用时可能需要更严格的验证。
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"The path {path} does not exist.")
    return os.path.abspath(path)

def get_xyz_filenames(directory: str) -> set:
    """
    获取目录中的所有.xyz文件名（不包含扩展名）。
    """
    directory = validate_path(directory)
    return {os.path.splitext(f)[0] for f in os.listdir(directory) if f.endswith('.gjf')}

def filter_csv_chunk(chunk: pd.DataFrame, xyz_filenames: set) -> pd.DataFrame:
    """
    过滤出Hash列值在xyz文件名中的数据。
    """
    return chunk[chunk['Hash'].isin(xyz_filenames)]

def filter_csv_by_xyz_files(csv_file_path: str, xyz_filenames: set, output_file_path: str, chunksize: int = DEFAULT_CHUNK_SIZE):
    """
    分块读取CSV文件，过滤出包含在xyz文件名中的数据，并保存到新的CSV文件。
    """
    try:
        filtered_chunks = []
        with ProcessPoolExecutor() as executor:
            # 使用并行处理来提高性能
            future_to_chunk = {executor.submit(filter_csv_chunk, chunk, xyz_filenames): chunk for chunk in pd.read_csv(csv_file_path, chunksize=chunksize)}
            for future in future_to_chunk:
                filtered_chunks.append(future.result())

        # 合并所有过滤后的块并保存到新的CSV文件
        filtered_df = pd.concat(filtered_chunks)
        filtered_df.to_csv(output_file_path, index=False)
        print(f"Filtered data saved to {output_file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

def main(source_directory: str, csv_file_path: str, output_file_path: str):
    """
    主函数，程序的入口点。
    """
    xyz_filenames = get_xyz_filenames(source_directory)
    filter_csv_by_xyz_files(csv_file_path, xyz_filenames, output_file_path)

if __name__ == '__main__':
    source_directory = validate_path(r'C:\Users\szou5\Desktop\work\gjf\unmatched')
    csv_file_path = validate_path(r'C:\Users\szou5\Desktop\work\results24.csv')
    output_file_path = validate_path(r'C:\Users\szou5\Desktop\work\results25.csv')
    main(source_directory, csv_file_path, output_file_path)