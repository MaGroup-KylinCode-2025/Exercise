#include "molecular.h"
<<<<<<< HEAD
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>

// 元素周期表数据结构
std::unordered_map<std::string, std::size_t> periodic_table = {
    {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
    {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},
    {"S", 16}, {"Cl", 17}, {"Ar", 18}
};

Molecular::Molecular(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // 读取电荷和自旋多重度
    file >> charge_ >> spin_;
    
    // 读取原子信息
    std::string element;
    double x, y, z;
    while (file >> element >> x >> y >> z) {
        Atom atom;
        atom.type = element;
        // 查找元素在周期表中的原子序数
        auto it = periodic_table.find(element);
        if (it != periodic_table.end()) {
            atom.Z = it->second;
        } else {
            atom.Z = 0; // 未知元素
        }
        atom.x = x;
        atom.y = y;
        atom.z = z;
        atoms_.push_back(atom);
    }
}

Molecular::Molecular(const Molecular &other) {
    atoms_ = other.atoms_;
    charge_ = other.charge_;
    spin_ = other.spin_;
    basis_ = other.basis_;
}

Molecular::Molecular(Molecular &&other) {
    atoms_ = std::move(other.atoms_);
    charge_ = other.charge_;
    spin_ = other.spin_;
    basis_ = std::move(other.basis_);
    
    other.charge_ = 0;
    other.spin_ = 0;
}

// Getters实现
const std::vector<Atom> Molecular::get_atoms() const {
    return atoms_;
}

const int Molecular::get_charge() const {
    return charge_;
}

const int Molecular::get_spin() const {
    return spin_;
}

const std::string Molecular::get_basis() const {
    return basis_;
}

// Setters实现
void Molecular::set_atoms(const std::vector<Atom> &atoms) {
    atoms_ = atoms;
}

void Molecular::set_charge(const int charge) {
    charge_ = charge;
}

void Molecular::set_spin(const int spin) {
    spin_ = spin;
}

void Molecular::set_basis(const std::string &basis) {
    basis_ = basis;
}

// 重载输出运算符
std::ostream &operator<<(std::ostream &os, const Molecular &molecule) {
    os << "=== Molecular Information ===" << std::endl;
    os << "Charge: " << molecule.charge_ << std::endl;
    os << "Spin multiplicity: " << molecule.spin_ << std::endl;
    os << "Basis set: " << molecule.basis_ << std::endl;
    os << "Number of atoms: " << molecule.atoms_.size() << std::endl;
    os << "Atoms:" << std::endl;
    
    for (const auto &atom : molecule.atoms_) {
        os << "  " << atom.type << " (Z=" << atom.Z << ") "
           << "at (" << atom.x << ", " << atom.y << ", " << atom.z << ")" << std::endl;
    }
    
    return os;
}
=======

Molecular::Molecular(const std::string &filename) {
  /*
  TODO: read filename and parse the file
  */
}
Molecular::Molecular(const Molecular &other) {
  /*
  TODO: copy constructor
  */
}
Molecular::Molecular(Molecular &&other) {
  /*
  TODO: move constructor
  */
}

>>>>>>> upstream/AAA-Source-tyb-branch
