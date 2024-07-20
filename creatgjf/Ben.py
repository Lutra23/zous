from rdkit import Chem
from rdkit.Chem import AllChem
import shutil

def generate_aromatic_structure(ring_smiles, n_rings):
    """
    根据共享边缘连接环生成多环芳烃结构。
    
    :param ring_smiles: 单个环的SMILES表示。
    :param n_rings: 环的数量。
    :return: 生成的多环芳烃结构作为RDKit分子对象。
    """
    # 验证n_rings是否为正整数
    if not isinstance(n_rings, int) or n_rings <= 0:
        raise ValueError("n_rings必须是正整数。")
    
    # 从SMILES创建第一个环并验证
    first_ring = Chem.MolFromSmiles(ring_smiles)
    if first_ring is None:
        raise ValueError(f"无法识别的环SMILES：{ring_smiles}")
    
    # 使用第一个环初始化多环分子
    aromatic_structure = Chem.RWMol(first_ring)
    
    # 当前分子的总原子数
    total_atoms = len(aromatic_structure.GetAtoms())
    
    # 准备一个模板环以提高效率
    template_ring = Chem.RWMol(first_ring)
    
    for i in range(1, n_rings):
        # 使用模板环并进行必要的调整
        new_ring = Chem.RWMol(template_ring)
        
        # 在插入前验证新环
        if new_ring is None:
            raise ValueError(f"错误创建环：为环{i+1}的模板修改失败。")
        
        aromatic_structure.InsertMol(new_ring)
        
        new_atoms = len(new_ring.GetAtoms())
        
        # 根据原子数调整连接索引
        connection_indices = [total_atoms - new_atoms + 1, total_atoms] if i != 1 else [1, 6]
        
        # 在新环和前一环之间添加单键连接
        aromatic_structure.AddBond(connection_indices[0], connection_indices[1], Chem.BondType.SINGLE)
        
        total_atoms += new_atoms
    
    # 确保所有环保持其芳香性并清理分子
    Chem.SanitizeMol(aromatic_structure)
    
    # 添加氢原子并嵌入分子
    aromatic_structure = Chem.AddHs(aromatic_structure)
    AllChem.EmbedMolecule(aromatic_structure)
    
    return aromatic_structure

def save_as_xyz(mol, filename):
    """
    将RDKit分子对象保存为XYZ文件。
    
    :param mol: RDKit分子对象。
    :param filename: 输出XYZ文件的名称。
    """
    # 验证输入是否为有效的RDKit分子对象
    if not isinstance(mol, Chem.Mol):
        raise ValueError("输入不是一个有效的RDKit分子对象。")
    
    # 验证文件名是否为字符串且以'.xyz'结尾
    if not isinstance(filename, str) or not filename.endswith('.xyz'):
        raise ValueError("无效的文件名。它必须是一个以'.xyz'结尾的字符串")
    
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    
    with open(filename, 'w') as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write("由RDKit生成\n")
        for atom in atoms:
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

if __name__ == '__main__':
    ring_smiles = 'c1ccccc1'
    output_dir = r'C:\Users\zous\Desktop\AIben\benzof_'  # 指定输出目录

    for n_rings in range(2, 7):
        try:
            multi_ring_mol = generate_aromatic_structure(ring_smiles, n_rings)
            #print(f"生成了{n_rings}环的多环芳烃结构。")
            
            xyzfile_name = f'{n_rings}.xyz'
            save_as_xyz(multi_ring_mol, xyzfile_name)
            #print(f"XYZ文件已保存为{xyzfile_name}")

            # 移动XYZ文件到指定目录
            shutil.move(xyzfile_name, output_dir + xyzfile_name)
            print(f"XYZ文件已移动到{output_dir}{xyzfile_name}")
        except ValueError as e:
            print(f"错误：{e}，无法生成{n_rings}环的结构")