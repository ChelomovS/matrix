#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include <type_traits>

#include "buffer.hpp"
#include "double_compare.hpp" 

namespace Matrix {

template <typename ElemT> 
class matrix_t final : private matrix_buffer_t<ElemT> {
    using matrix_buffer_t<ElemT>::rows_;
    using matrix_buffer_t<ElemT>::cols_;
    using matrix_buffer_t<ElemT>::buf_;
//===================================================================================================
    private:
    struct ProxyRow {
        ElemT *row;

        const ElemT& operator[](int n) const {
            return row[n];
        }

        ElemT& operator[](int n) {
            return row[n];
        }
    };
//===================================================================================================
    public:
    int get_rows() const {
        return rows_; 
    }

    int get_cols() const {
        return cols_;
    }

    ElemT multiply_diag() const {
        ElemT result = 1;

        for (auto i = 0; i < rows_; ++i) {
            result *= buf_[i * cols_ + i];
        }

        return result;
    }

    ElemT trace() const {
        ElemT result = 0;

        for (auto i = 0; i < rows_; ++i) {
            result += buf_[i * cols_ + i];
        }

        return result;
    }

    ProxyRow operator[](int n) const {
        return ProxyRow{buf_ + n * cols_};
    }

    ElemT get_det_by_gauss_algorithm() const {
        int rows = get_rows();
        int cols = get_cols();

        matrix_t<double> double_matrix{*this}; 

        double det = 1.0; 

        for (int i = 0; i < rows; ++i) {
            int k = i;
            for (int j = i + 1; j < rows; ++j) {
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

            for (int j = i + 1; j < rows; ++j) {
                double coeff = double_matrix[j][i] / double_matrix[i][i];
                for (int c = i; c < cols; ++c) {
                    double_matrix[j][c] -= coeff * double_matrix[i][c];
                }
            }
        }

        return std::is_floating_point_v<ElemT> ? det : std::round(det);
    }
//===================================================================================================
    public:
    void swap_rows(int first_row_number, int second_row_number) {
        if (first_row_number == second_row_number) 
            return;

        ElemT* first_row  = buf_ + first_row_number  * cols_;
        ElemT* second_row = buf_ + second_row_number * cols_;

        for (auto i = 0; i < cols_; ++i) {
            std::swap(first_row[i], second_row[i]);
        }
    }

    matrix_t& negate() & {
        for (auto i = 0; i < rows_; ++i) {
            for (auto j = 0; j < cols_; ++j) {
                buf_[i * cols_ + j] *= -1;
            }
        }

        return *this;
    }

    matrix_t& transpose() & {
        matrix_t transposed{cols_, rows_};
        for (auto i = 0; i < rows_; ++i) {
            for (auto j = 0; j < cols_; ++j) {
                transposed.buf_[j * rows_ + i] = buf_[i * cols_ + j];
            }
        }

        std::swap(buf_, transposed.buf_);
        std::swap(rows_, cols_);
        
        return *this;
    }
//===================================================================================================
    public:
    matrix_t(int rows, int cols): matrix_buffer_t<ElemT>{rows, cols} {}

    matrix_t(int rows, int cols, ElemT value): matrix_buffer_t<ElemT>{rows, cols} {
        for (auto i = 0; i < rows_ * cols_; ++i) {
            Constructor(buf_ + i, value);
        } 
    }

    template <typename Iterator>
    matrix_t(int rows, int cols, Iterator start, Iterator end): matrix_buffer_t<ElemT>{rows, cols} {
        int number_of_elem = 0;

        ElemT elem{};

        for (auto i = start; i < end; ++i) {
            elem = static_cast<ElemT>(*i);
            Constructor(buf_ + number_of_elem, elem);
            ++number_of_elem;
        }
    }

    matrix_t(const matrix_t& other): matrix_buffer_t<ElemT>{other.get_cols(), 
                                                            other.get_rows()} {
        for (auto i = 0; i < rows_ * cols_; ++i) {
            Constructor(buf_ + i, other.buf_[i]);
        }
    }

    template <typename AnotherElemT> explicit
    matrix_t(const matrix_t<AnotherElemT>& other): matrix_buffer_t<ElemT>{other.get_cols(), 
                                                                          other.get_rows()} {
        for (auto i = 0; i < rows_; ++i) {
            for (auto j = 0; j < cols_; ++j) {
                buf_[i * rows_ + j] = static_cast<ElemT>(other[i][j]);
            }
        }
    }

    matrix_t& operator=(const matrix_t& other) { 
        matrix_t<ElemT> tmp{other};
        std::swap(*this, tmp);
        return *this;
    }

    matrix_t(matrix_t&& other) noexcept   = default;
    matrix_t& operator=(matrix_t&& other) = default;

    ~matrix_t() {
        Destructor(buf_, buf_ + rows_ * cols_);
    }
};
//=================================================================================================
template <typename ElemT>
bool operator==(const matrix_t<ElemT>& lhs, const matrix_t<ElemT>& rhs) {
    if ((lhs.get_cols() != rhs.get_cols()) || (lhs.get_rows() != rhs.get_rows()))
        return false;

    for (auto i = 0; i < lhs.get_rows(); ++i) {
        for (auto j = 0; j < lhs.get_cols(); ++j) {
            if (!Compare::is_equal(lhs[i][j], rhs[i][j]))
                return false;
        }
    }

    return true;
}

template <typename ElemT> 
std::istream& operator>>(std::istream& inp_stream, matrix_t<ElemT>& matrix) {
    int rows = matrix.get_rows();
    int cols = matrix.get_cols();

    for (auto i = 0; i < rows; ++i) {
        for (auto j = 0; j < cols; ++j) {
            inp_stream >> matrix[i][j];
        }
    }

    return inp_stream;
}

template <typename ElemT> 
std::ostream& operator<<(std::ostream& out_stream, const matrix_t<ElemT>& matrix) {
    int rows = matrix.get_rows();
    int cols = matrix.get_cols();
    
    for (auto i = 0; i < rows; ++i) {
        for (auto j = 0; j < cols; ++j) {
            out_stream << matrix[i][j] << " ";
        }

        out_stream << std::endl;
    }

    out_stream << std::endl;

    return out_stream;
}

template <typename ElemT>
bool operator!=(const matrix_t<ElemT>& lhs, const matrix_t<ElemT>& rhs) {
    return !(lhs == rhs);
}

} // namespace Matrix 

#endif // MATRIX_HPP
