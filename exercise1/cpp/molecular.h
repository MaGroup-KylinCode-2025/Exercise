#ifndef MOLECULAR_H
#define MOLECULAR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>

// 元素周期表
static const std::unordered_map<std::string, std::size_t> PeriodicTable{
    {"X", 0}, {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4},
    {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9},
    // ... (完整的周期表如文档所示)
};

// 单位转换常量
constexpr double Angstrom2AU{1.0 / 0.5291772083};
constexpr double AU2Angstrom{0.5291772083};

// 原子结构体
struct Atom {
    std::string symbol; // 元素符号
    std::size_t Z;      // 原子序数
    double x, y, z;     // 坐标
};

class Molecular {
public:
    // 构造函数
    Molecular() = default;
    explicit Molecular(const std::string &filename);
    Molecular(const Molecular &other) noexcept;
    Molecular(Molecular &&other) noexcept;
    
    // Getter函数
    inline int get_charge() const noexcept;
    inline int get_spin() const noexcept;
    inline const std::vector<Atom> &get_atoms() const noexcept;
    inline const std::string &get_basis() const noexcept;
    
    // Setter函数
    inline void set_charge(const int charge) noexcept;
    inline void set_spin(const int spin) noexcept;
    inline void set_atoms(const std::vector<Atom> &atoms) noexcept;
    inline void set_basis(const std::string &basis) noexcept;
    
    // 运算符重载
    friend std::ostream &operator<<(std::ostream &os, const Molecular &molecule);

private:
    std::vector<Atom> atoms_;           // 原子列表
    int charge_;                        // 电荷
    int spin_;                          // 自旋多重度
    std::string basis_;                 // 基组名称
    std::array<std::size_t, 3> nelec_;  // [alpha电子数, beta电子数, 总电子数]
};

#endif // MOLECULAR_H
