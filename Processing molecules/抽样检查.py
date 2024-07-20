import pandas as pd
import random
from rdkit import Chem
from rdkit.Chem import Draw
import logging
import os

# 设置日志记录
logging.basicConfig(level=logging.INFO)

# 确保RDKit的绘图功能可用
if not os.environ.get('_RDKIT年报图_', ''):
    import matplotlib
    matplotlib.use('Agg')


def smiles_to_png(smiles, output_file, image_size=(300, 300)):
    """
    将 SMILES 表示转换为 PNG 图片
    
    参数：
        smiles (str): SMILES 表示的分子结构
        output_file (str): 输出 PNG 图片的文件名
        image_size (tuple): 输出图像的宽度和高度，默认为(300, 300)
        
    """

    try:
        # 将 SMILES 转换为分子对象
        molecule = Chem.MolFromSmiles(smiles)
        
        if molecule is None:
            logging.error("无法从SMILES字符串创建分子对象: {}".format(smiles))
            raise ValueError("无效的SMILES字符串")

        # 设置图像大小
        Draw.MolToFile(molecule, output_file, size=image_size)
        logging.info("成功将SMILES转换为PNG图像: {}".format(output_file))

    except Exception as e:
        logging.error("在转换过程中遇到错误: {}".format(e))
        raise
# 假设CSV文件有两列：'ID'和'SMILES'
def sample_and_draw_from_csv(csv_path, sample_size=10, output_dir='output_images'):
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取CSV文件
    df = pd.read_csv(csv_path)
    
    # 确保SMILES列存在
    if 'SMILES' not in df.columns:
        raise ValueError("CSV文件中未找到'SMILES'列")
    
    # 随机抽样
    sampled_smiles = df['SMILES'].sample(sample_size).tolist()
    
    for idx, smiles in enumerate(sampled_smiles, start=1):
        # 使用之前定义的函数转换并保存图片
        output_file = os.path.join(output_dir, f'sample_{idx}.png')
        try:
            smiles_to_png(smiles, output_file)
            print(f"成功转换并保存了样本 {idx} 至 {output_file}")
        except Exception as e:
            print(f"转换样本 {idx} 时出错: {e}")

# 调用函数，传入你的CSV路径和抽样大小
csv_path = 'filtered_molecules.csv'  # 替换为你的CSV文件路径
sample_and_draw_from_csv(csv_path, sample_size=20)