#include "molecular.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>

// 有参构造函数
Molecular::Molecular(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // 读取原子数
    int atom_count;
    infile >> atom_count;
    atoms_.reserve(atom_count);
    
    // 跳过第一行剩余部分
    std::string line;
    std::getline(infile, line);
    
    // 读取第二行：charge 和 spin
    int charge, spin;
    infile >> charge >> spin;
    charge_ = charge;
    spin_ = spin;
    std::getline(infile, line); // 跳过第二行剩余部分
    
    // 计算总电子数
    int total_nelec = -charge;
    
    // 读取原子信息
    for (int i = 0; i < atom_count; ++i) {
        std::string symbol;
        double x, y, z;
        if (!(infile >> symbol >> x >> y >> z)) {
            throw std::runtime_error("Invalid atom line in file: " + filename);
        }
        
        // 格式化元素符号
        if (symbol.size() == 1) {
            symbol[0] = toupper(symbol[0]);
        } else if (symbol.size() == 2) {
            symbol[0] = toupper(symbol[0]);
            symbol[1] = tolower(symbol[1]);
        } else {
            throw std::runtime_error("Invalid atom symbol: " + symbol);
        }
        
        auto Z = PeriodicTable.at(symbol);
        atoms_.push_back({symbol, Z, x * Angstrom2AU, y * Angstrom2AU, z * Angstrom2AU});
        total_nelec += Z;
    }
    
    // 初始化电子数数组 [alpha, beta, total]
    nelec_[2] = total_nelec;  // 总电子数
    nelec_[0] = (total_nelec + spin) / 2;  // alpha电子数
    nelec_[1] = total_nelec - nelec_[0];   // beta电子数
    
    // 检查自旋设置是否合理
    if ((total_nelec + spin) % 2 != 0) {
        throw std::runtime_error("Invalid spin multiplicity for given charge and atoms");
    }
    
    basis_ = ""; // 默认为空
}

// 拷贝构造函数
Molecular::Molecular(const Molecular &other) noexcept
    : atoms_(other.atoms_), charge_(other.charge_), spin_(other.spin_),
      basis_(other.basis_), nelec_(other.nelec_) {}

// 移动构造函数  
Molecular::Molecular(Molecular &&other) noexcept
    : atoms_(std::move(other.atoms_)), charge_(other.charge_),
      spin_(other.spin_), basis_(std::move(other.basis_)), nelec_(other.nelec_) {}

// Getter函数实现
inline int Molecular::get_charge() const noexcept { return charge_; }
inline int Molecular::get_spin() const noexcept { return spin_; }
inline const std::vector<Atom> &Molecular::get_atoms() const noexcept { return atoms_; }
inline const std::string &Molecular::get_basis() const noexcept { return basis_; }

// Setter函数实现
inline void Molecular::set_charge(const int charge) noexcept { charge_ = charge; }
inline void Molecular::set_spin(const int spin) noexcept { spin_ = spin; }
inline void Molecular::set_atoms(const std::vector<Atom> &atoms) noexcept { atoms_ = atoms; }
inline void Molecular::set_basis(const std::string &basis) noexcept { basis_ = basis; }

// 输出运算符重载
std::ostream &operator<<(std::ostream &os, const Molecular &molecule) {
    os << "Molecular object\n"
       << " Charge: " << molecule.get_charge() << "\n"
       << " Spin: " << molecule.get_spin() << "\n"
       << " Basis: " << molecule.get_basis() << "\n"
       << " Atoms: " << molecule.get_atoms().size() << "\n";
    for (const auto &atom : molecule.get_atoms()) {
        os << "  " << atom.symbol << " (" << atom.Z << "): " 
           << atom.x << " " << atom.y << " " << atom.z << "\n";
    }
    return os;
}
