from pathlib import Path


PeriodicTable = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
}

Angstrom2AU = 1.8897259886


class Molecular:
    def __init__(self, file_name: str):
        self.file_name = Path(file_name)
        self.atoms_ = []
        self.charge_ = 0
        self.spin_ = 0
        self.basis_ = "cc-pVDZ"
        self.nelec_ = (0, 0, 0)  # (total, alpha, beta)
        self._initialize()

    def _initialize(self):
        if not self.file_name.exists():
            raise FileNotFoundError(f"Cannot open file: {self.file_name}")

        with open(self.file_name, "r") as f:
            atom_count = int(f.readline().split("#", 1)[0].strip())
            charge_spin_line = f.readline().split("#", 1)[0].strip()
            charge, spin = map(int, charge_spin_line.split())
            self.charge_ = charge
            self.spin_ = spin

            total_nelec = -charge
            for _ in range(atom_count):
                line = f.readline().split("#", 1)[0].strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) != 4:
                    raise ValueError(f"Invalid atom line in file: {self.file_name}")
                symbol, x, y, z = parts
                symbol = symbol.capitalize()
                if symbol not in PeriodicTable:
                    raise ValueError(f"Unknown atom symbol: {symbol}")

                Z = PeriodicTable[symbol]
                x, y, z = (
                    float(x) * Angstrom2AU,
                    float(y) * Angstrom2AU,
                    float(z) * Angstrom2AU,
                )
                self.atoms_.append({"symbol": symbol, "Z": Z, "x": x, "y": y, "z": z})
                total_nelec += Z

            alpha = (spin + total_nelec) // 2
            beta = total_nelec - alpha
            self.nelec_ = (total_nelec, alpha, beta)

            self._check_nelec(total_nelec, alpha, beta, spin)

    def _check_nelec(self, total, alpha, beta, spin):
        if total != alpha + beta:
            raise ValueError(
                f"Invalid number of electrons: total={total}, alpha={alpha}, beta={beta}"
            )
        if spin != alpha - beta:
            raise ValueError(f"Invalid spin: spin={spin}, alpha={alpha}, beta={beta}")

    # print the information of the molecule
    def print(self):
        print("Molecular object")
        print(f"  Charge: {self.charge_}")
        print(f"  Spin:   {self.spin_}")
        print(f"  Basis:  {self.basis_}")
        print(f"  Atoms:  {len(self.atoms_)}")
        for atom in self.atoms_:
            print(
                f"    {atom['symbol']} ({atom['Z']}): {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}"
            )

    # Getters
    def get_atoms(self):
        return self.atoms_

    def get_charge(self):
        return self.charge_

    def get_spin(self):
        return self.spin_

    def get_basis(self):
        return self.basis_

    # Setters
    def set_atoms(self, atoms):
        self.atoms_ = atoms

    def set_charge(self, charge):
        self.charge_ = charge

    def set_spin(self, spin):
        self.spin_ = spin

    def set_basis(self, basis):
        self.basis_ = basis


if __name__ == "__main__":
    mol = Molecular("water.xyz")
    mol.print()
