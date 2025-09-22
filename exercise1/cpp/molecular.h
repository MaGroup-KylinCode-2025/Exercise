#ifndef MOLECULAR_H
#define MOLECULAR_H

#include <iostream>
#include <source_location>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

static const std::unordered_map<std::string, std::size_t> PeriodicTable{
    {"X", 0},    {"H", 1},    {"He", 2},   {"Li", 3},   {"Be", 4},
    {"B", 5},    {"C", 6},    {"N", 7},    {"O", 8},    {"F", 9},
    {"Ne", 10},  {"Na", 11},  {"Mg", 12},  {"Al", 13},  {"Si", 14},
    {"P", 15},   {"S", 16},   {"Cl", 17},  {"Ar", 18},  {"K", 19},
    {"Ca", 20},  {"Sc", 21},  {"Ti", 22},  {"V", 23},   {"Cr", 24},
    {"Mn", 25},  {"Fe", 26},  {"Co", 27},  {"Ni", 28},  {"Cu", 29},
    {"Zn", 30},  {"Ga", 31},  {"Ge", 32},  {"As", 33},  {"Se", 34},
    {"Br", 35},  {"Kr", 36},  {"Rb", 37},  {"Sr", 38},  {"Y", 39},
    {"Zr", 40},  {"Nb", 41},  {"Mo", 42},  {"Tc", 43},  {"Ru", 44},
    {"Rh", 45},  {"Pd", 46},  {"Ag", 47},  {"Cd", 48},  {"In", 49},
    {"Sn", 50},  {"Sb", 51},  {"Te", 52},  {"I", 53},   {"Xe", 54},
    {"Cs", 55},  {"Ba", 56},  {"La", 57},  {"Ce", 58},  {"Pr", 59},
    {"Nd", 60},  {"Pm", 61},  {"Sm", 62},  {"Eu", 63},  {"Gd", 64},
    {"Tb", 65},  {"Dy", 66},  {"Ho", 67},  {"Er", 68},  {"Tm", 69},
    {"Yb", 70},  {"Lu", 71},  {"Hf", 72},  {"Ta", 73},  {"W", 74},
    {"Re", 75},  {"Os", 76},  {"Ir", 77},  {"Pt", 78},  {"Au", 79},
    {"Hg", 80},  {"Tl", 81},  {"Pb", 82},  {"Bi", 83},  {"Po", 84},
    {"At", 85},  {"Rn", 86},  {"Fr", 87},  {"Ra", 88},  {"Ac", 89},
    {"Th", 90},  {"Pa", 91},  {"U", 92},   {"Np", 93},  {"Pu", 94},
    {"Am", 95},  {"Cm", 96},  {"Bk", 97},  {"Cf", 98},  {"Es", 99},
    {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104},
    {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109},
    {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114},
    {"Lv", 116}, {"Ts", 117}, {"Og", 118}};

constexpr double Angstrom2AU{1.0 / 0.5291772083};
constexpr double AU2Angstrom{0.5291772083};

struct Atom {
  std::string symbol; // element
  std::size_t Z;      // atomic number
  double x, y, z;     // coordinates
};

class Molecular {
public:
  // Default constructor
  Molecular() = default;
  // Constructor with filename
  Molecular(const std::string &filename);
  // Copy constructor
  Molecular(const Molecular &other) noexcept;
  // Move constructor
  Molecular(Molecular &&other) noexcept;
  // Destructor
  ~Molecular() = default;

  // Getters
  inline std::span<const Atom> get_atoms() const noexcept;
  inline int get_charge() const noexcept;
  inline int get_spin() const noexcept;
  inline const std::string &get_basis() const noexcept;

  // Setters
  inline void set_atoms(const std::vector<Atom> &atoms) noexcept;
  inline void set_charge(const int charge) noexcept;
  inline void set_spin(const int spin) noexcept;
  inline void set_basis(const std::string &basis) noexcept;

  // Help function to check spin and nelec_;
  void check_nelec_(
      const std::array<std::size_t, 3> &nelec, std::size_t spin,
      const std::source_location &loc = std::source_location::current());

  // Overload << operator
  friend std::ostream &operator<<(std::ostream &os, const Molecular &molecule);

private:
  std::vector<Atom> atoms_;          // atoms in the molecule
  int charge_;                       // charge of the molecule
  int spin_;                         // spin multiplicity of the molecule
  std::string basis_;                // basis set name
  std::array<std::size_t, 3> nelec_; // number of alpha/beta electrons
}; // class Molecular
#endif // MOLECULAR_H