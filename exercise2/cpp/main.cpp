#include "timer.hpp"
#include <cblas.h>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>

template <typename T>
std::vector<T> matrix_multiply(const std::vector<T> &a, const std::vector<T> &b,
                               const int ld) {
  Timer<> timer;
  auto row = a.size() / ld;
  auto col = b.size() / ld;
  std::vector<T> c(row * col, T(0));
  for (auto i = 0; i < row; ++i) {
    for (auto j = 0; j < col; ++j) {
      for (auto k = 0; k < ld; ++k) {
        // c_ij += a_ik * b_kj
        c[i * col + j] += a[i * ld + k] * b[k * col + j];
      }
    }
  }
  return c;
}

template <typename T>
std::vector<T> matrix_multiply_blas(const std::vector<T> &a,
                                    const std::vector<T> &b, const int ld) {
  Timer<> timer;
  auto row = a.size() / ld;
  auto col = b.size() / ld;
  std::vector<T> c(row * col, T(0));
  if constexpr (std::is_same_v<T, float>) {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, 1.0,
                a.data(), ld, b.data(), col, 0.0, c.data(), col);
  } else if constexpr (std::is_same_v<T, double>) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, 1.0,
                a.data(), ld, b.data(), col, 0.0, c.data(), col);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    const float alpha[2] = {1.0f, 0.0f};
    const float beta[2] = {0.0f, 0.0f};
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, 
                alpha, a.data(), ld, b.data(), col, beta, c.data(), col);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    const double alpha[2] = {1.0, 0.0};
    const double beta[2] = {0.0, 0.0};
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, 
                alpha, a.data(), ld, b.data(), col, beta, c.data(), col);
  } else {
    throw std::runtime_error("Unsupported type");
  }
  return c;
}

template <typename T>
void print_matrix(const std::vector<T> &m, const int row, const int col) {
  if (row * col != m.size()) {
    throw std::runtime_error("Invalid matrix size");
  }
  for (auto i = 0; i < row; ++i) {
    for (auto j = 0; j < col; ++j) {
      T val = m[i * col + j];
      if constexpr (std::is_same_v<T, float>) {
        std::cout << std::fixed << std::setprecision(2) << val << " ";
      } else if constexpr (std::is_same_v<T, double>) {
        std::cout << std::fixed << std::setprecision(2) << val << " ";
      } else if constexpr (std::is_same_v<T, std::complex<float>> || 
                          std::is_same_v<T, std::complex<double>>) {
        std::cout << std::fixed << std::setprecision(2) 
                  << "(" << val.real() << "," << val.imag() << ") ";
      }
    }
    std::cout << "\n";
  }
}

int main() {
  std::cout << "=== 矩阵乘法性能测试 ===\n" << std::endl;

  // 测试1: 小矩阵验证正确性
  std::cout << "测试1: 小矩阵 (2x2) 验证正确性" << std::endl;
  std::vector<float> a = {1.0f, 2.0f, 3.0f, 4.0f};
  std::vector<float> b = {5.0f, 6.0f, 7.0f, 8.0f};
  
  std::cout << "矩阵 A:" << std::endl;
  print_matrix(a, 2, 2);
  std::cout << "矩阵 B:" << std::endl;
  print_matrix(b, 2, 2);
  
  std::cout << "\n朴素实现结果:" << std::endl;
  auto c1 = matrix_multiply(a, b, 2);
  print_matrix(c1, 2, 2);
  
  std::cout << "\nBLAS实现结果:" << std::endl;
  auto c2 = matrix_multiply_blas(a, b, 2);
  print_matrix(c2, 2, 2);
  
  // 测试2: 大矩阵性能对比
  std::cout << "\n\n测试2: 大矩阵性能对比" << std::endl;
  const int N = 500;  // 矩阵大小
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis(0.0f, 1.0f);
  
  std::vector<float> A(N * N);
  std::vector<float> B(N * N);
  for (int i = 0; i < N * N; ++i) {
    A[i] = dis(gen);
    B[i] = dis(gen);
  }
  
  std::cout << "矩阵大小: " << N << "x" << N << std::endl;
  
  std::cout << "\n朴素实现 (Naive Implementation):" << std::endl;
  {
    auto C1 = matrix_multiply(A, B, N);
  }
  
  std::cout << "\nBLAS实现 (OpenBLAS):" << std::endl;
  {
    auto C2 = matrix_multiply_blas(A, B, N);
  }
  
  // 测试3: 不同数据类型
  std::cout << "\n\n测试3: 不同数据类型 (double)" << std::endl;
  std::vector<double> Ad(N * N);
  std::vector<double> Bd(N * N);
  std::uniform_real_distribution<double> dis_d(0.0, 1.0);
  for (int i = 0; i < N * N; ++i) {
    Ad[i] = dis_d(gen);
    Bd[i] = dis_d(gen);
  }
  
  std::cout << "矩阵大小: " << N << "x" << N << " (double)" << std::endl;
  std::cout << "BLAS实现:" << std::endl;
  {
    auto Cd = matrix_multiply_blas(Ad, Bd, N);
  }
  
  std::cout << "\n=== 测试完成 ===" << std::endl;
  
  return 0;
}
