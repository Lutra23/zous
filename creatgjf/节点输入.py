import os
import glob

def merge_gjf_files(output_file, input_folder):
    try:
        # 搜索输入文件夹中的所有 .gjf 文件
        input_files = glob.glob(os.path.join(input_folder, '*.gjf'))
        
        if not input_files:
            print(f"文件夹 {input_folder} 中没有找到 .gjf 文件。")
            return
        
        with open(output_file, 'w') as outfile:
            for i, fname in enumerate(input_files):
                with open(fname, 'r') as infile:
                    outfile.write(infile.read())
                if i < len(input_files) - 1:
                    outfile.write('\n--Link1--\n')
        
        print(f"所有文件已成功合并到 {output_file}")
    except Exception as e:
        print(f"合并文件时发生错误: {e}")

if __name__ == "__main__":
    # 输出文件名
    output_file = "output.gjf"
    
    # 输入文件夹路径
    input_folder = r"C:\Users\zous\Desktop\NHC"
    
    # 调用合并函数
    merge_gjf_files(output_file, input_folder)
