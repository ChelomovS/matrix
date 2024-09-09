#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>

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
    long long get_rows() const {
        return rows_; 
    }

    long long get_cols() const {
        return cols_;
    }

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

    T det() const {
        
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

    //---------------------------------------------------------------------------------
    public:
    Matrix(long long rows, long long cols): Buffer<T>{rows, cols} {}

    Matrix(const Matrix<T> &other): Buffer<T>{other.get_rows(), other.get_cols()} {
        long long number_of_elem = 0;
        for (auto i = 0; i < rows_; ++i) {
            for (auto j = 0; j < cols_; ++j) {
                Constructor(buf_ + number_of_elem, static_cast<T>(other[i][j]));
            }
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
    //---------------------------------------------------------------------------------
    private:
    void make_good_colomn(long long position_in_col) const { //REVIEW -  move to another place
        for (auto i = position_in_col + 1; i < rows_; ++i) {
            double pivot_elem = buf_[position_in_col][position_in_col];
            double mult_coef  = buf_[i][j] / pivot_elem;
            if (is_equal(mult_coef, 0)) 
                continue;
            for (auto j = position_in_col; j < cols_; ++j) {
                //TODO - !!!
            }
        }
    }

    long long max_abs_elem_in_col(Matrix<double>& matrix, long long col) const {
        long long position       = 0;
        double    max_elem       = 0;
        long long number_of_rows = matrix.get_rows();

        for (auto i = 0; i < number_of_rows; ++i) {
            double next_elem = std::fabs(matrix[i][col]);
            if (next_elem > max_elem) {
                max_elem = next_elem;
                position = i;
            }
        }
        return position;
    }
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