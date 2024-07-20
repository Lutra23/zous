from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import logging
import shutil

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def generate_polybenzene(n_rings):
    """
    生成指定数量苯环的聚苯乙烯结构。
    """
    if not isinstance(n_rings, int) or n_rings <= 0:
        raise ValueError("n_rings 必须是正整数")
    
    benzene = Chem.MolFromSmiles('c1ccccc1')
    if benzene is None:
        logging.error("无法从SMILES字符串生成苯的分子结构")
        raise Exception("生成苯的分子结构失败")
    
    polybenzene = Chem.RWMol(benzene)

    try:
        for i in range(1, n_rings):
            new_benzene = Chem.MolFromSmiles('c1ccccc1')
            if new_benzene is None:
                logging.error(f"无法从SMILES字符串生成苯环 {i+1} 的分子结构")
                continue
            
            polybenzene.InsertMol(new_benzene)
            
            # 连接新的苯环到现有结构
            polybenzene.AddBond((i - 1) * 6 + 1, i * 6 + 4, Chem.BondType.SINGLE)

        Chem.SanitizeMol(polybenzene)
        logging.info("分子结构成功标准化")

        polybenzene = Chem.AddHs(polybenzene)
        if AllChem.EmbedMolecule(polybenzene) != 0:
            logging.error("无法嵌入分子")
            raise Exception("分子嵌入失败")
        
        logging.info("分子结构成功添加氢原子并嵌入")
    except Exception as e:
        logging.error(f"生成聚苯乙烯结构时出现错误: {str(e)}")
        raise

    return polybenzene

def convert_to_smiles(mol):
    """
    将RDKit分子对象转换为SMILES字符串。
    """
    return Chem.MolToSmiles(mol, kekuleSmiles=False)

def save_as_xyz(mol, filename):
    """
    将RDKit分子对象保存为XYZ文件。
    """
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    
    try:
        with open(filename, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\n")
            f.write("Generated by RDKit\n")
            for atom in atoms:
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n")
            logging.info(f"成功保存XYZ文件至 {filename}")
    except Exception as e:
        logging.error(f"保存XYZ文件时出现错误: {str(e)}")

if __name__ == '__main__':
    output_dir = r'C:\Users\zous\Desktop\AIben\polybenzene_'  # 指定输出目录

    for n_rings in range(2, 7):
        try:
            polybenzene_mol = generate_polybenzene(n_rings)

            smiles = convert_to_smiles(polybenzene_mol)
            print(f'SMILES: {smiles}')

            xyzfile_name = f'{n_rings}.xyz'
            save_as_xyz(polybenzene_mol, xyzfile_name)
            print(f'XYZ文件保存为 {xyzfile_name}')

            shutil.move(xyzfile_name, output_dir + xyzfile_name)
            print(f"XYZ文件已移动到{output_dir}{xyzfile_name}")
        except ValueError as e:
            print(f"错误：{e}，无法生成{n_rings}环的结构")