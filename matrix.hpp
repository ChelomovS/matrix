#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>

#include <cmath>
#include "buffer.hpp"
#include "double_compare.hpp" 

namespace Matrix {

template <typename T = double> 
class Matrix final : private Buffer<T> {
    using Buffer<T>::rows_;
    using Buffer<T>::cols_;
    using Buffer<T>::buf_;
    //---------------------------------------------------------------------------------
    private:
    struct ProxyRow {
        T *row;
        const T& operator[](int n) const {
            return row[n];
        }

        T& operator[](int n) {
            return row[n];
        }
    };
    //---------------------------------------------------------------------------------
    public:    
    ProxyRow operator[](int n) const {
        if (n < 0 || n > rows_) {
            std::cerr << "Invalid param!" << std::endl;
            std::abort();
        }
        return ProxyRow{buf_ + n * cols_};
    }

    T multiply_diag() const {
        if (rows_ != cols_) {
            std::cerr << "Matrix is not square!" << std::endl;
            std::abort();
        }
        T result = 1;

        for (auto i = 0; i < rows_; ++i) {
            result *= buf_[i][i];
        }
        return result;
    }
    long long get_rows() const {
        return rows_; 
    }

    long long get_cols() const {
        return cols_;
    }

    Matrix& negative(const Matrix& matrix) {
        Matrix tmp_matrix{matrix};
        rows = matrix.get_rows();
        cols = matrix.get_cols();
        for (auto i = 0; i < rows; ++i) {
            for (auto j = 0; j < cols; ++j) {
                tmp_matrix *= -1;
            }
        }
        return tmp_matrix;
    }

    T trace() const {
        if (rows_ != cols_) {
            std::cerr << "Matrix is not square!" << std::endl;
            std::abort();
        }
        T result = 0;
        for (auto i = 0; i < rows_; ++i) {
            result += buf_[i][i];
        }
        return result;
    }

    T det(const Matrix& matrix) {
        Matrix tmp_matrix{matrix};
        return det_by_gauss_algorithm(tmp_matrix);
    }

    // Complexity: O(N^3)
    T det_by_gauss_algorithm(Matrix& matrix) {
        long long rows = matrix.get_rows();
        long long cols = matrix.get_cols();
        if (rows != cols) {
            std::cerr << "Matrix is not square" << std::endl;
            std::abort();
        }

        T det = 0;

        for (auto i = 0; i < rows; ++i) { 
            long long k = i; 
            // Find max elem in column 
            for (auto j = i + 1; j < rows; ++j) {
                if (abs(matrix[j][i])) > abs(matrix[k][j])
                    k = j;
            }

            if (is_equal(matrix[k][i], 0)) {
                // If a row in a matrix consists of zero elements, 
                // then the determinant is zero
                det = 0;
                return det;
            }
            if (i != k) {
                std::swap(matrix[i], matrix[k]);
                // If you change two different rows in a matrix,
                // the sign of the determinant changes to the opposite
                det *= -1; 
            }

            det *= matrix[i][i]
            for (auto j = i + 1; j < rows; ++j) {
                matrix[i][j] /= matrix[i][i];
            }

            for (auto j = 0; j < rows; ++j) {
                if ((j != i) && (!is_equal(matrix[j][i], 0))) {
                    for (auto k = i + 1; k < n; k++) {
                        matrix[j][k] -= matrix[i][k] * matrix[j][i];
                    }
                }
            }
        }
        return det;
    }
    //---------------------------------------------------------------------------------
    public:
    Matrix(long long rows, long long cols): Buffer<T>{rows, cols} {}
    Matrix(long long rows, long long cols, T val =T{}): Buffer<T>(rows, cols) {
        long long number_of_elem = 0;
        for (auto i = 0; i < rows; ++i) {
            for (auto j = 0; j < cols; ++j) {
                Сonstructor(buf_ + number_of_elem, val)ж
            }
        }
    }
    Matrix(const Matrix<T> &other): Buffer<T>{other.get_rows(), other.get_cols()} {
        long long number_of_elem = 0;
        for (auto i = 0; i < rows_; ++i) {
            for (auto j = 0; j < cols_; ++j) {
                Constructor(buf_ + number_of_elem, static_cast<T>(other[i][j]));
                ++number_of_elem;
            }
        }
    }
    template <typename It>
    Matrix(long long rows, long long cols, It start, It end): Buffer<T>{rows, cols} {
        long long number_of_elem = 0;
        T elem;
        for (It i = start; i < end; ++i) {
            elem = static_cast<T>(*i);
            Constructor(buf_ + number_of_elem, elem);
            ++number_of_elem;
        }
    }
    Matrix& operator=(const Matrix &other) {
        Matrix<T> tmp{other};
        std::swap(*this, tmp);
        return *this;
    }

    Matrix(Matrix &&other) noexcept  = default;
    Matrix& operator=(Matrix &&other) = default;

    ~Matrix() = default; 
};

//-------------------------------------------------------------------------------------
template <typename T>
bool operator==(const Matrix<T> lhs, const Matrix<T> rhs) {
    if ((lhs.get_cols() != rhs.get_cols()) || (lhs.get_rows() != rhs.get_rows()))
        return false;
    for (auto i = 0; i < lhs.get_rows(); ++i) {
        for (auto j = 0; j < lhs.get_cols(); ++j){
            if(!is_equal(lhs[i][j], rhs[i][j]))
                return false;
        }
    }
    return true;
}

template <typename T> 
std::istream& operator>>(std::istream& inp_stream, Matrix<T>& matrix) {
    long long rows = matrix.get_rows();
    long long cols = matrix.get_cols();
    for (auto i = 0; i < rows; ++i) {
        for (auto j = 0; j < cols; ++j) {
            inp_stream >> matrix[i][j];
        }
    }
    return inp_stream;
}

template <typename T> 
std::ostream& operator<<(std::ostream& out_stream, Matrix<T>& matrix) {
    long long rows = matrix.get_rows();
    long long cols = matrix.get_cols();
    for (auto i = 0; i < rows; ++i) {
        for (auto j = 0; j < cols; ++j){
            out_stream << matrix[i][j] << " ";
        }
        out_stream << std::endl;
    }
    out_stream << std::endl;
    return out_stream;
}

template <typename T>
bool operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    return !(lhs == rhs);
}

} // namespace Matrix 

#endif // MATRIX_HPP