import re
import shutil
from pathlib import Path

class GJFChecker:
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.charge = None
        self.spin = None
        self.atoms = []

    def read_gjf(self):
        with self.filepath.open('r', encoding='utf-8') as file:
            lines = file.readlines()

        charge_spin_line_found = False
        for line in lines:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            if not charge_spin_line_found:
                if re.match(r'^\s*-?\d+\s+\d+\s*$', line):
                    self.charge, self.spin = map(int, line.split())
                    charge_spin_line_found = True
            else:
                atom_data = line.split()
                if len(atom_data) >= 4:
                    self.atoms.append(atom_data)

    def validate(self):
        if self.charge is None or self.spin is None:
            raise ValueError("未找到电荷和自旋多重度。")
        electron_count = 0
        for atom in self.atoms:
            element = atom[0]
            if element not in atomic_numbers:
                raise ValueError(f"未知元素: {element}")
            electron_count += atomic_numbers[element]

        electron_count -= self.charge
        expected_spin = 2 * ((electron_count % 2) == 1) + 1

        if self.spin != expected_spin:
            raise ValueError(f"不合理的自旋多重度: {self.spin}。根据电荷 {self.charge} 和电子数 {electron_count}，预期为 {expected_spin}。")

    def check(self):
        self.read_gjf()
        self.validate()

# 完整的元素原子序数列表
atomic_numbers = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
    'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
    'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
    'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84,
    'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94,
    'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
    'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
    'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}

def process_gjf_files(directory):
    dir_path = Path(directory)
    if not dir_path.is_dir():
        raise ValueError("指定路径不是文件夹。")

    reasonable_dir = dir_path / 'reasonable'
    unreasonable_dir = dir_path / 'unreasonable'
    reasonable_dir.mkdir(exist_ok=True)
    unreasonable_dir.mkdir(exist_ok=True)

    for gjf_file in dir_path.glob('*.gjf'):
        checker = GJFChecker(gjf_file)
        try:
            checker.check()
            shutil.move(str(gjf_file), reasonable_dir / gjf_file.name)
            print(f"{gjf_file.name}: 合理")
        except ValueError as e:
            shutil.move(str(gjf_file), unreasonable_dir / gjf_file.name)
            print(f"{gjf_file.name}: 不合理 - {e}")

if __name__ == "__main__":
    directory = r'C:\Users\szou5\Desktop\nn ns\XYZ\unreasonable'  # 修改为你的GJF文件所在文件夹路径
    try:
        process_gjf_files(directory)
        print("处理完成。")
    except ValueError as e:
        print(e)

