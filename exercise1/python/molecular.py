class Molecular:
    def __init__(self, file_name):
<<<<<<< HEAD
        self.charge_ = 0
        self.spin_ = 0
        self.basis_ = None
        self.atoms_ = []
        self.initial_lize(file_name)

    def initial_lize(self, file_name):
        """从XYZ文件读取原子信息"""
        try:
            with open(file_name, 'r') as f:
                lines = f.readlines()
                # 跳过前两行（原子数和注释行）
                for line in lines[2:]:
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        atom_symbol = parts[0]
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        self.atoms_.append({
                            'symbol': atom_symbol,
                            'x': x,
                            'y': y,
                            'z': z
                        })
        except FileNotFoundError:
            print(f"错误：文件 {file_name} 未找到")
        except Exception as e:
            print(f"读取文件时出错：{e}")

    def print(self):
        """打印分子信息"""
        print("=" * 50)
        print("分子信息：")
        print(f"电荷: {self.charge_}")
        print(f"自旋: {self.spin_}")
        print(f"基组: {self.basis_}")
        print("原子列表：")
        for i, atom in enumerate(self.atoms_):
            print(f"  原子 {i+1}: {atom['symbol']} "
                  f"({atom['x']:.6f}, {atom['y']:.6f}, {atom['z']:.6f})")
        print("=" * 50)

    # Getter 方法
=======
        # TODO:
        self.charge_ = 0
        self.spin_ = 0
        self.basis_ = None
        self.atoms_ = None
        self.initial_lize()

    def initial_lize(self):
        # TODO:
        pass

    # print the information of the molecule
    def print(self):
        pass

    # Getter
>>>>>>> upstream/AAA-Source-tyb-branch
    def get_atoms(self):
        return self.atoms_

    def get_charge(self):
        return self.charge_

    def get_spin(self):
        return self.spin_

    def get_basis(self):
        return self.basis_

<<<<<<< HEAD
    # Setter 方法
=======
    # Setter
>>>>>>> upstream/AAA-Source-tyb-branch
    def set_atoms(self, atoms):
        self.atoms_ = atoms

    def set_charge(self, charge):
<<<<<<< HEAD
        self.charge_ = charge  # 修正拼写错误
=======
        self.charge_ = charge
>>>>>>> upstream/AAA-Source-tyb-branch

    def set_spin(self, spin):
        self.spin_ = spin

    def set_basis(self, basis):
        self.basis_ = basis


if __name__ == "__main__":
<<<<<<< HEAD
    # 创建示例XYZ文件（如果不存在）
    xyz_content = """3
水分子示例
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000"""
    
    with open("h2o.xyz", "w") as f:
        f.write(xyz_content)
    
    # 测试分子类
    mol = Molecular("h2o.xyz")
    mol.set_charge(0)      # 设置电荷
    mol.set_spin(1)        # 设置自旋
    mol.set_basis("6-31G") # 设置基组
    mol.print()
    
    # 测试getter方法
    print(f"获取电荷: {mol.get_charge()}")
    print(f"获取自旋: {mol.get_spin()}")
    print(f"获取基组: {mol.get_basis()}")
    print(f"获取原子数: {len(mol.get_atoms())}")
=======
    mol = Molecular("h2o.xyz")
    mol.print()
>>>>>>> upstream/AAA-Source-tyb-branch
