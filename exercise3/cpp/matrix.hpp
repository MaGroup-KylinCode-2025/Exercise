#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cblas.h>
#include <iostream>
#include <lapacke.h>
#include <list>
#include <omp.h>
#include <stdexcept>
#include <complex>
#include <vector>
#include <cstring>

template <typename T> class Matrix {
public:
  Matrix() : rows_(0), cols_(0), data_(nullptr), order_(CblasRowMajor) {}

  Matrix(int rows, int cols, CBLAS_ORDER order = CblasRowMajor)
      : rows_(rows), cols_(cols), data_(new T[rows * cols]), order_(order) {
    std::memset(data_, 0, sizeof(T) * rows * cols);
  }
  
  Matrix(int rows, int cols, T *data)
      : rows_(rows), cols_(cols), data_(data),
        order_(CBLAS_ORDER::CblasRowMajor) {
    if (data == nullptr) {
      throw std::invalid_argument("data is nullptr");
    }
  }
  
  Matrix(int rows, int cols, std::initializer_list<T> list)
      : rows_(rows), cols_(cols), data_(new T[rows * cols]),
        order_(CblasRowMajor) {
    if (list.size() != rows * cols) {
      throw std::invalid_argument("list size not match");
    }
    std::copy(list.begin(), list.end(), data_);
  }

  ~Matrix() {
    if (data_ != nullptr) {
      delete[] data_;
    }
  }

  // 拷贝构造函数
  Matrix(const Matrix &other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    order_ = other.order_;
    data_ = new T[rows_ * cols_];
    std::memcpy(data_, other.data_, sizeof(T) * rows_ * cols_);
  }
  
  // 移动构造函数
  Matrix(Matrix &&other) noexcept {
    rows_ = other.rows_;
    cols_ = other.cols_;
    order_ = other.order_;
    data_ = other.data_;
    other.data_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }

  const int rows() const { return rows_; }
  const int cols() const { return cols_; }
  const T *data() const { return data_; }
  T *data() { return data_; }

  // operator() - 元素访问
  T &operator()(int row, int col) {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
      throw std::out_of_range("Matrix subscript out of range");
    }
    if (order_ == CblasRowMajor) {
      return data_[row * cols_ + col];
    } else {
      return data_[col * rows_ + row];
    }
  }
  
  const T &operator()(int row, int col) const {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
      throw std::out_of_range("Matrix subscript out of range");
    }
    if (order_ == CblasRowMajor) {
      return data_[row * cols_ + col];
    } else {
      return data_[col * rows_ + row];
    }
  }

  // operator + - 矩阵加法
  Matrix<T> operator+(const Matrix<T> &other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::invalid_argument("Matrix size not match");
    }
    Matrix<T> result(rows_, cols_);
    int n = rows_ * cols_;
    
    if constexpr (std::is_same_v<T, float>) {
      cblas_scopy(n, data_, 1, result.data_, 1);
      cblas_saxpy(n, 1.0f, other.data_, 1, result.data_, 1);
    } else if constexpr (std::is_same_v<T, double>) {
      cblas_dcopy(n, data_, 1, result.data_, 1);
      cblas_daxpy(n, 1.0, other.data_, 1, result.data_, 1);
    } else {
      // 对于复数类型，使用手动循环
      for (int i = 0; i < n; ++i) {
        result.data_[i] = data_[i] + other.data_[i];
      }
    }
    return result;
  }
  
  // operator - - 矩阵减法
  Matrix<T> operator-(const Matrix<T> &other) const {
    if (rows_ != other.rows() || cols_ != other.cols()) {
      throw std::invalid_argument("Matrix size not match");
    }
    Matrix<T> result(rows_, cols_);
    int n = rows_ * cols_;
    
    if constexpr (std::is_same_v<T, float>) {
      cblas_scopy(n, data_, 1, result.data_, 1);
      cblas_saxpy(n, -1.0f, other.data_, 1, result.data_, 1);
    } else if constexpr (std::is_same_v<T, double>) {
      cblas_dcopy(n, data_, 1, result.data_, 1);
      cblas_daxpy(n, -1.0, other.data_, 1, result.data_, 1);
    } else {
      for (int i = 0; i < n; ++i) {
        result.data_[i] = data_[i] - other.data_[i];
      }
    }
    return result;
  }

  // operator * - 矩阵乘法
  Matrix<T> operator*(const Matrix<T> &other) const {
    if (cols_ != other.rows()) {
      throw std::invalid_argument("Matrix size not match");
    }
    Matrix<T> result(rows_, other.cols());
    
    if constexpr (std::is_same_v<T, float>) {
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                  rows_, other.cols(), cols_, 
                  1.0f, data_, cols_, other.data_, other.cols(),
                  0.0f, result.data_, result.cols());
    } else if constexpr (std::is_same_v<T, double>) {
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  rows_, other.cols(), cols_,
                  1.0, data_, cols_, other.data_, other.cols(),
                  0.0, result.data_, result.cols());
    } else if constexpr (std::is_same_v<T, std::complex<float>>) {
      std::complex<float> alpha(1.0f, 0.0f), beta(0.0f, 0.0f);
      cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  rows_, other.cols(), cols_,
                  &alpha, data_, cols_, other.data_, other.cols(),
                  &beta, result.data_, result.cols());
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
      std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  rows_, other.cols(), cols_,
                  &alpha, data_, cols_, other.data_, other.cols(),
                  &beta, result.data_, result.cols());
    }
    return result;
  }
  
  // operator * - 标量乘法
  Matrix<T> operator*(T scalar) const {
    Matrix<T> result(rows_, cols_);
    int n = rows_ * cols_;
    
    if constexpr (std::is_same_v<T, float>) {
      cblas_scopy(n, data_, 1, result.data_, 1);
      cblas_sscal(n, scalar, result.data_, 1);
    } else if constexpr (std::is_same_v<T, double>) {
      cblas_dcopy(n, data_, 1, result.data_, 1);
      cblas_dscal(n, scalar, result.data_, 1);
    } else {
      for (int i = 0; i < n; ++i) {
        result.data_[i] = data_[i] * scalar;
      }
    }
    return result;
  }
  
  // operator *= - 标量乘法（原地）
  Matrix<T> &operator*=(T scalar) {
    int n = rows_ * cols_;
    if constexpr (std::is_same_v<T, float>) {
      cblas_sscal(n, scalar, data_, 1);
    } else if constexpr (std::is_same_v<T, double>) {
      cblas_dscal(n, scalar, data_, 1);
    } else {
      for (int i = 0; i < n; ++i) {
        data_[i] *= scalar;
      }
    }
    return *this;
  }

  // operator << - 打印矩阵
  friend std::ostream &operator<<(std::ostream &os, const Matrix<T> &mat) {
    for (int i = 0; i < mat.rows_; ++i) {
      for (int j = 0; j < mat.cols_; ++j) {
        os << mat(i, j) << " ";
      }
      os << "\n";
    }
    return os;
  }

private:
  int rows_, cols_;
  T *data_;
  CBLAS_ORDER order_;
};

// 类型信息辅助结构
template <typename T> struct MatrixTypeInfo {
  using RealType = T;
};
template <> struct MatrixTypeInfo<std::complex<float>> {
  using RealType = float;
};
template <> struct MatrixTypeInfo<std::complex<double>> {
  using RealType = double;
};

// eigh - 标准本征值问题 HC = eC
template <typename T>
std::pair<std::vector<typename MatrixTypeInfo<T>::RealType>, Matrix<T>>
eigh(const Matrix<T> &H) {
  using RealType = typename MatrixTypeInfo<T>::RealType;
  if (H.rows() != H.cols()) {
    throw std::invalid_argument("Matrix must be square");
  }
  int n = H.rows();
  std::vector<RealType> w(n);
  Matrix<T> eigvec(H);

  int info = 0;
  if constexpr (std::is_same_v<T, float>) {
    info = LAPACKE_ssyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n, w.data());
  } else if constexpr (std::is_same_v<T, double>) {
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n, w.data());
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    info = LAPACKE_cheevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n, w.data());
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    info = LAPACKE_zheevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n, w.data());
  } else {
    throw std::runtime_error("Unsupported type for eigh");
  }

  if (info != 0) {
    throw std::runtime_error("LAPACK diagonalization failed");
  }
  return {w, eigvec};
}

// eigh - 广义本征值问题 HC = eSC
template <typename T>
std::pair<std::vector<typename MatrixTypeInfo<T>::RealType>, Matrix<T>>
eigh(const Matrix<T> &H, const Matrix<T> &S) {
  using RealType = typename MatrixTypeInfo<T>::RealType;
  if (H.rows() != H.cols() || S.rows() != S.cols()) {
    throw std::invalid_argument("Matrices must be square");
  }
  if (H.rows() != S.rows()) {
    throw std::invalid_argument("Matrix sizes must match");
  }
  
  int n = H.rows();
  std::vector<RealType> w(n);
  Matrix<T> H_copy(H);
  Matrix<T> S_copy(S);

  int info = 0;
  if constexpr (std::is_same_v<T, float>) {
    info = LAPACKE_ssygvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', n, 
                          H_copy.data(), n, S_copy.data(), n, w.data());
  } else if constexpr (std::is_same_v<T, double>) {
    info = LAPACKE_dsygvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', n,
                          H_copy.data(), n, S_copy.data(), n, w.data());
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    info = LAPACKE_chegvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', n,
                          H_copy.data(), n, S_copy.data(), n, w.data());
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    info = LAPACKE_zhegvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', n,
                          H_copy.data(), n, S_copy.data(), n, w.data());
  } else {
    throw std::runtime_error("Unsupported type for eigh");
  }

  if (info != 0) {
    throw std::runtime_error("LAPACK generalized eigenvalue problem failed");
  }
  return {w, H_copy};
}

#endif // MATRIX_HPP
