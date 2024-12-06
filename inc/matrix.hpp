#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <type_traits>

#include "buffer.hpp"
#include "double_compare.hpp" 

namespace Matrix {

template <typename ElemT> 
class matrix_t final : private matrix_buffer_t<ElemT> {
    using matrix_buffer_t<ElemT>::buf_ ;
    using matrix_buffer_t<ElemT>::rows_;
    using matrix_buffer_t<ElemT>::cols_;
    using matrix_buffer_t<ElemT>::used_;
//===================================================================================================
    private:
    struct ProxyRow {
        ElemT *row;

        const ElemT& operator[](size_t n) const {
            return row[n];
        }

        ElemT& operator[](size_t n) {
            return row[n];
        }
    };
//===================================================================================================
    public:
    size_t get_used() const noexcept {
        return used_;
    }

    size_t get_rows() const noexcept {
        return rows_; 
    }

    size_t get_cols() const noexcept {
        return cols_;
    }

    ElemT multiply_diag() const {
        ElemT result = 1;

        for (size_t i = 0; i < rows_; ++i) {
            result *= buf_[i * cols_ + i];
        }

        return result;
    }

    ElemT trace() const {
        ElemT result = 0;

        for (size_t i = 0; i < rows_; ++i) {
            result += buf_[i * cols_ + i];
        }

        return result;
    }

    ProxyRow operator[](size_t n) const {
        return ProxyRow{buf_ + n * cols_};
    }

    ElemT get_det_by_gauss_algorithm() const {
        matrix_t<double> double_matrix{*this}; 

        double det = 1.0; 

        for (size_t i = 0; i < rows_; ++i) {
            size_t k = i;
            for (size_t j = i + 1; j < rows_; ++j) {
                if (std::abs(double_matrix[j][i]) > std::abs(double_matrix[k][i])) {
                    k = j;
                }
            }

            if (Compare::is_equal(double_matrix[k][i], 0.0)) {
                return 0.0;
            }

            if (i != k) {
                double_matrix.swap_rows(i, k);
                det *= -1.0;
            }

            det *= double_matrix[i][i];

            for (size_t j = i + 1; j < rows_; ++j) {
                double coeff = double_matrix[j][i] / double_matrix[i][i];
                for (size_t c = i; c < cols_; ++c) {
                    double_matrix[j][c] -= coeff * double_matrix[i][c];
                }
            }
        }

        return std::is_floating_point_v<ElemT> ? det : std::round(det);
    }
//===================================================================================================
    public:
    void swap_rows(size_t first_row_number, size_t second_row_number) {
        if (first_row_number == second_row_number) 
            return;

        ElemT* first_row  = buf_ + first_row_number  * cols_;
        ElemT* second_row = buf_ + second_row_number * cols_;

        for (size_t i = 0; i < cols_; ++i) {
            std::swap(first_row[i], second_row[i]);
        }
    }

    matrix_t& negate() & {
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                buf_[i * cols_ + j] *= -1;
            }
        }

        return *this;
    }

    matrix_t& transpose() & {
        matrix_t transposed{cols_, rows_};
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                transposed.buf_[j * rows_ + i] = buf_[i * cols_ + j];
            }
        }

        std::swap(buf_, transposed.buf_);
        std::swap(rows_, cols_);
        
        return *this;
    }
//===================================================================================================
    public:
    explicit matrix_t(size_t rows, size_t cols): matrix_buffer_t<ElemT>{rows, cols} {}

    explicit matrix_t(size_t rows, size_t cols, const ElemT& value): matrix_buffer_t<ElemT>{rows, cols} {
        while (used_ < rows_ * cols_) {
            Constructor(buf_ + used_, value);
            ++used_;
        } 
    }

    template <typename Iterator>
    matrix_t(size_t rows, size_t cols, Iterator start, Iterator end): matrix_t<ElemT>{rows, cols} {
        ElemT elem{};

        for (auto iter = start; iter != end; ++iter) {
            elem = static_cast<ElemT>(*iter);
            Constructor(buf_ + used_, elem);
            ++used_;
        }
    }

    template <typename AnotherElemT> explicit
    matrix_t(const matrix_t<AnotherElemT>& other): matrix_buffer_t<ElemT>{other.get_rows(), other.get_cols()} {
        size_t current_row = 0;
        size_t current_col = 0;

        while (used_ < other.get_used()) {
            Constructor(buf_ + used_, static_cast<ElemT>(other[current_row][current_col]));
            ++used_;
            if (current_row == rows_ - 1) {
                current_row = 0;
                ++current_col;
            } else {
                ++current_row;
            }
        }
    }
    
    matrix_t(const matrix_t& other): matrix_buffer_t<ElemT>{other.get_rows(), other.get_cols()} {
        size_t current_row = 0;
        size_t current_col = 0;

        while (used_ < other.get_used()) {
            Constructor(buf_ + used_, other[current_row][current_col]);
            ++used_;
            if (current_row == rows_ - 1) {
                current_row = 0;
                ++current_col;
            } else {
                ++current_row;
            }
        }
    }

    matrix_t& operator=(const matrix_t& other) { 
        matrix_t<ElemT> tmp{other};
        std::swap(*this, tmp);
        return *this;
    }

    matrix_t(matrix_t&& other) noexcept            = default;
    matrix_t& operator=(matrix_t&& other) noexcept = default;

    ~matrix_t() = default;
};
//=================================================================================================
template <typename ElemT>
bool operator==(const matrix_t<ElemT>& lhs, const matrix_t<ElemT>& rhs) {
    if ((lhs.get_cols() != rhs.get_cols()) || (lhs.get_rows() != rhs.get_rows()))
        return false;

    for (size_t i = 0; i < lhs.get_rows(); ++i) {
        for (size_t j = 0; j < lhs.get_cols(); ++j) {
            if (!Compare::is_equal(lhs[i][j], rhs[i][j]))
                return false;
        }
    }

    return true;
}

template <typename ElemT>
bool operator!=(const matrix_t<ElemT>& lhs, const matrix_t<ElemT>& rhs) {
    return !(lhs == rhs);
}

template <typename ElemT> 
std::istream& operator>>(std::istream& inp_stream, matrix_t<ElemT>& matrix) {
    size_t rows = matrix.get_rows();
    size_t cols = matrix.get_cols();

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            inp_stream >> matrix[i][j];
        }
    }

    return inp_stream;
}

template <typename ElemT> 
std::ostream& operator<<(std::ostream& out_stream, const matrix_t<ElemT>& matrix) {
    size_t rows = matrix.get_rows();
    size_t cols = matrix.get_cols();
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            out_stream << matrix[i][j] << " ";
        }

        out_stream << std::endl;
    }

    out_stream << std::endl;

    return out_stream;
}

} // namespace Matrix 

#endif // MATRIX_HPP
