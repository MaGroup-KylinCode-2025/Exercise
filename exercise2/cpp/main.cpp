#include "timer.hpp"
#include <cblas.h>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

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
                a.data(), ld, b.data(), ld, 0.0, c.data(), ld);
  } else if constexpr (std::is_same_v<T, double>) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, 1.0,
                a.data(), ld, b.data(), ld, 0.0, c.data(), ld);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    std::complex<float> alpha(1.0f, 0.0f);
    std::complex<float> beta(0.0f, 0.0f);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, &alpha,
                a.data(), ld, b.data(), ld, &beta, c.data(), ld);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    std::complex<double> alpha(1.0, 0.0);
    std::complex<double> beta(0.0, 0.0);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, &alpha,
                a.data(), ld, b.data(), ld, &beta, c.data(), ld);
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
      }
      // TODO: Add other types
    }
    std::cout << "\n";
  }
}

int main() {
  std::vector<float> a = {
      1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
  };
  auto b = a;
  int ld = 3;
  auto c = matrix_multiply<float>(a, b, ld);
  std::cout << "Matrix multiplication (naive):\n";
  print_matrix(c, ld, ld);

  auto c_blas = matrix_multiply_blas<float>(a, b, ld);
  std::cout << "Matrix multiplication (BLAS):\n";
  print_matrix(c_blas, ld, ld);

  for (auto i = 0; i < ld * ld; ++i) {
    if (std::abs(c[i] - c_blas[i]) > 1e-5) {
      std::cout << "BLAS result is not correct\n";
      break;
    }
  }

  return 0;
}