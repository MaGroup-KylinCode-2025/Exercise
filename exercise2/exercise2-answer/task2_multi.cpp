#include <iostream>
#include "task2_timer.hpp" //typename T = std::chrono::nanoseconds
#include <vector>
#include <complex>
#include <stdexcept>
#include <chrono>
#include <cblas.h> // OpenBLAS

// Matrix multiplication using nested loops
template <typename T>
std::vector<T> matrix_multiply(const std::vector<T> &a, const std::vector<T> &b, const int ld) {
    auto row = a.size() / ld;
    auto col = b.size() / ld;
    std::vector<T> c(row * col, T(0));
    for (auto i = 0; i < row; ++i) {
        for (auto j = 0; j < col; ++j) {
            for (auto k = 0; k < ld; ++k) {
                c[i * col + j] += a[i * ld + k] * b[k * col + j];
            }
        }
    }
    return c;
}

// Matrix multiplication using BLAS
template <typename T>
std::vector<T> matrix_multiply_blas(const std::vector<T> &a, const std::vector<T> &b, const int ld) {
    // Timer<> timer; // Measure time
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
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, &T(1.0),
                    reinterpret_cast<const float*>(a.data()), ld, reinterpret_cast<const float*>(b.data()), ld, &T(0.0), reinterpret_cast<float*>(c.data()), ld);
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, &T(1.0),
                    reinterpret_cast<const double*>(a.data()), ld, reinterpret_cast<const double*>(b.data()), ld, &T(0.0), reinterpret_cast<double*>(c.data()), ld);
    } else {
        throw std::runtime_error("Unsupported type");
    }
    return c;
}

int main() {
    // Example usage
    int ld = 3; // Leading dimension
    std::vector<double> a = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> b = {9, 8, 7, 6, 5, 4, 3, 2, 1};

    // Measure time for nested loop multiplication
    {
        Timer<> timer; // Measure time
        auto result = matrix_multiply(a, b, ld);
        std::cout << "Nested loop result: ";
        for (const auto &val : result) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    // Measure time for BLAS multiplication
    {
        Timer<> timer; // Measure time
        auto result = matrix_multiply_blas(a, b, ld);
        std::cout << "BLAS result: ";
        for (const auto &val : result) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}