from rdkit import Chem

def generate_benzenes_with_rings(num_rings):
    benzenes = []
    for _ in range(num_rings):
        if num_rings == 1:
            benzene_smiles = "c1ccccc1"
        elif num_rings >= 2:
            benzene_smiles = "c1" + "ccccc" * (num_rings - 1) + "c1"
        benzene = Chem.MolFromSmiles(benzene_smiles)
        benzenes.append(benzene)
    return benzenes

# 生成 3 个具有不同环数量的并苯分子的 SMILES 表示法
num_rings_list = [1, 2, 3]
for num_rings in num_rings_list:
    benzenes = generate_benzenes_with_rings(num_rings)
    print(f"Benzenes with {num_rings} rings:")
    for i, benzene in enumerate(benzenes):
        if benzene:
            print(f"Benzene {i+1} SMILES:", Chem.MolToSmiles(benzene))
        else:
            print(f"Error generating benzene {i+1} with {num_rings} rings.")
    print()
