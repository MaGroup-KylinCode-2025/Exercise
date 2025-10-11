#include "matrix.hpp"
#include <iostream>
#include <iomanip>

// Huckel理论计算丁二烯（butadiene）的π电子结构
int main() {
  std::cout << "=== Huckel Molecular Orbital Theory ===" << std::endl;
  std::cout << "Example: Butadiene (CH2=CH-CH=CH2)" << std::endl << std::endl;

  // 丁二烯有4个π电子（4个p轨道）
  // Huckel哈密顿矩阵：
  // - 对角元：α（原子能量，通常设为0）
  // - 相邻原子：β（共振积分，通常设为-1）
  // - 非相邻原子：0
  
  // 构造4x4的Huckel矩阵（α=0, β=-1）
  int n = 4;  // 4个碳原子
  Matrix<double> H(n, n, {
    0.0, -1.0,  0.0,  0.0,   // C1
   -1.0,  0.0, -1.0,  0.0,   // C2
    0.0, -1.0,  0.0, -1.0,   // C3
    0.0,  0.0, -1.0,  0.0    // C4
  });

  std::cout << "Hamiltonian Matrix H:" << std::endl;
  std::cout << H << std::endl;

  // 对角化得到本征值和本征向量
  auto [eigenvalues, eigenvectors] = eigh(H);

  std::cout << "Eigenvalues (Orbital Energies in units of β):" << std::endl;
  for (int i = 0; i < n; ++i) {
    std::cout << "E" << i+1 << " = " << std::setw(10) << std::fixed 
              << std::setprecision(6) << eigenvalues[i] << " β" << std::endl;
  }

  std::cout << "\nEigenvectors (Molecular Orbital Coefficients):" << std::endl;
  std::cout << "Each column is a molecular orbital" << std::endl;
  std::cout << eigenvectors << std::endl;

  // 计算π电子总能量
  // 丁二烯有4个π电子，填充在最低的2个轨道中（每个轨道2个电子）
  double total_energy = 2.0 * (eigenvalues[0] + eigenvalues[1]);
  std::cout << "Total π-electron energy = " << total_energy << " β" << std::endl;
  std::cout << "                        = " << total_energy * (-1.0) << " |β|" << std::endl;

  // 离域能：与两个孤立的乙烯分子相比
  double localized_energy = 4.0 * (-1.0);  // 2个乙烯：2*(2*(-1))
  double delocalization_energy = total_energy - localized_energy;
  std::cout << "\nDelocalization energy = " << delocalization_energy << " β" 
            << std::endl;
  std::cout << "                      = " << delocalization_energy * (-1.0) 
            << " |β|" << std::endl;

  // 另一个例子：苯（benzene）
  std::cout << "\n\n=== Example 2: Benzene (C6H6) ===" << std::endl;
  
  int n_benzene = 6;
  Matrix<double> H_benzene(n_benzene, n_benzene, {
    0.0, -1.0,  0.0,  0.0,  0.0, -1.0,
   -1.0,  0.0, -1.0,  0.0,  0.0,  0.0,
    0.0, -1.0,  0.0, -1.0,  0.0,  0.0,
    0.0,  0.0, -1.0,  0.0, -1.0,  0.0,
    0.0,  0.0,  0.0, -1.0,  0.0, -1.0,
   -1.0,  0.0,  0.0,  0.0, -1.0,  0.0
  });

  auto [e_benzene, c_benzene] = eigh(H_benzene);

  std::cout << "\nBenzene Orbital Energies:" << std::endl;
  for (int i = 0; i < n_benzene; ++i) {
    std::cout << "E" << i+1 << " = " << std::setw(10) << std::fixed 
              << std::setprecision(6) << e_benzene[i] << " β" << std::endl;
  }

  // 苯有6个π电子
  double benzene_energy = 2.0 * (e_benzene[0] + e_benzene[1] + e_benzene[2]);
  std::cout << "\nTotal π-electron energy = " << benzene_energy << " β" << std::endl;
  std::cout << "                        = " << benzene_energy * (-1.0) 
            << " |β|" << std::endl;

  // 芳香稳定化能
  double three_ethylene = 6.0 * (-1.0);
  double aromatic_stabilization = benzene_energy - three_ethylene;
  std::cout << "\nAromatic stabilization energy = " << aromatic_stabilization 
            << " β" << std::endl;
  std::cout << "                              = " 
            << aromatic_stabilization * (-1.0) << " |β|" << std::endl;

  return 0;
}
