#include "molecular.h"
#include <fstream>

Molecular::Molecular(const std::string &filename) {
  std::ifstream infile(filename);
  if (!infile) {
    throw std::runtime_error("Cannot open file: " + filename);
  }
  int atom_count;
  infile >> atom_count;
  std::string line;
  std::getline(infile, line);

  int charge, spin;
  infile >> charge >> spin;
  charge_ = charge;
  spin_ = spin;
  std::getline(infile, line);
  int total_nelec{-charge};
  atoms_.reserve(atom_count);
  for (int i = 0; i < atom_count; ++i) {
    std::string symbol;
    double x, y, z;
    if (!(infile >> symbol >> x >> y >> z)) {
      throw std::runtime_error("Invalid atom line in file: " + filename);
    }

    if (symbol.size() == 1) {
      symbol[0] = toupper(symbol[0]);
    } else if (symbol.size() == 2) {
      symbol[0] = toupper(symbol[0]);
      symbol[1] = tolower(symbol[1]);
    } else {
      throw std::runtime_error("Invalid atom symbol: " + symbol);
    }
    auto Z = PeriodicTable.at(symbol);
    atoms_.push_back(
        {symbol, Z, x * Angstrom2AU, y * Angstrom2AU, z * Angstrom2AU});
    total_nelec += Z;
  }

  basis_ = "cc-pVDZ";

  std::size_t alpha = (spin + total_nelec) / 2;
  std::size_t beta = total_nelec - alpha;
  nelec_ = {static_cast<std::size_t>(total_nelec), alpha, beta};

  check_nelec_(nelec_, spin_);
}

void Molecular::check_nelec_(const std::array<std::size_t, 3> &nelec,
                             std::size_t spin,
                             const std::source_location &loc) {
  if (nelec[0] != nelec[1] + nelec[2]) {
    throw std::runtime_error(
        "Invalid number of electrons at " + std::string(loc.file_name()) + ":" +
        std::to_string(loc.line()) + " (" + loc.function_name() + ")\n" +
        "  total=" + std::to_string(nelec[0]) + " alpha=" +
        std::to_string(nelec[1]) + " beta=" + std::to_string(nelec[2]));
  }

  if (spin != nelec[1] - nelec[2]) {
    throw std::runtime_error("Invalid spin at " + std::string(loc.file_name()) +
                             ":" + std::to_string(loc.line()) + " (" +
                             loc.function_name() + ")\n" +
                             "  spin=" + std::to_string(spin) +
                             " alpha=" + std::to_string(nelec[1]) +
                             " beta=" + std::to_string(nelec[2]));
  }
}

Molecular::Molecular(const Molecular &other) noexcept
    : atoms_(other.atoms_), charge_(other.charge_), spin_(other.spin_),
      basis_(other.basis_) {}

Molecular::Molecular(Molecular &&other) noexcept
    : atoms_(std::move(other.atoms_)), charge_(other.charge_),
      spin_(other.spin_), basis_(std::move(other.basis_)) {}

// Getters
inline std::span<const Atom> Molecular::get_atoms() const noexcept {
  return atoms_;
}
inline int Molecular::get_charge() const noexcept { return charge_; }
inline int Molecular::get_spin() const noexcept { return spin_; }
inline const std::string &Molecular::get_basis() const noexcept {
  return basis_;
}

// Setters
inline void Molecular::set_atoms(const std::vector<Atom> &atoms) noexcept {
  atoms_ = atoms;
}
inline void Molecular::set_charge(const int charge) noexcept {
  charge_ = charge;
}
inline void Molecular::set_spin(const int spin) noexcept { spin_ = spin; }
inline void Molecular::set_basis(const std::string &basis) noexcept {
  basis_ = basis;
}

// Overload << operator
std::ostream &operator<<(std::ostream &os, const Molecular &molecule) {
  os << "Molecular object\n"
     << "  Charge: " << molecule.get_charge() << "\n"
     << "  Spin:   " << molecule.get_spin() << "\n"
     << "  Basis:  " << molecule.get_basis() << "\n"
     << "  Atoms:  " << molecule.get_atoms().size() << "\n";
  for (const auto &atom : molecule.get_atoms()) {
    os << "    " << atom.symbol << " (" << atom.Z << "): " << atom.x << " "
       << atom.y << " " << atom.z << "\n";
  }
  return os;
}
