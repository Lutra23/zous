from rdkit import Chem
from rdkit.Chem import AllChem

def xyz_to_smiles(xyz_file_path):
    """
    将XYZ格式的三维坐标文件转换为SMILES字符串
    参数:
    - xyz_file_path: XYZ文件的路径

    返回:
    - smiles: 生成的SMILES字符串
    """
    # 使用RDKit从XYZ文件读取分子
    # 注意：RDKit原生不支持直接从.xyz文件读取，因此这里需要用其他方式处理
    # 这个例子假设XYZ文件格式是标准的，即第一行是原子数量，第二行是注释，后面是原子坐标
    with open(xyz_file_path, 'r') as file:
        lines = file.readlines()
    
    # 从文件内容构建一个分子
    mol = Chem.RWMol()
    for line in lines[2:]:  # 跳过前两行
        parts = line.split()
        if len(parts) < 4:
            continue  # 如果不是标准的坐标行，则跳过
        atom_symbol, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
        atom = Chem.Atom(atom_symbol)
        idx = mol.AddAtom(atom)
    
    # 添加虚拟键，以便RDKit能够识别分子结构
    for idx in range(mol.GetNumAtoms()):
        if idx > 0:
            mol.AddBond(idx-1, idx, Chem.BondType.SINGLE)
    
    # 生成三维坐标
    AllChem.EmbedMolecule(mol)
    
    # 生成二维坐标，以便于在ChemDraw中更好地呈现
    AllChem.Compute2DCoords(mol)

    # 生成SMILES字符串
    smiles = Chem.MolToSmiles(mol)
    
    return smiles

# 这里替换为你的XYZ文件路径
xyz_file_path = 'path/to/your/molecule.xyz'
smiles = xyz_to_smiles(xyz_file_path)
print("Generated SMILES:", smiles)