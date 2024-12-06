#ifndef BUFFER_HPP
#define BUFFER_HPP

#include <algorithm>
#include <utility>

namespace Matrix {
//===================================================================================================
template <typename ElemT> void Constructor(ElemT* ptr, const ElemT& value) {
    new (ptr) ElemT{value};
}

template <typename ElemT> void Destructor(ElemT* ptr) noexcept {
    ptr->~ElemT();
}

template <typename Iterator> void Destructor(Iterator begin, Iterator end) noexcept {
    while (begin != end) {
        Destructor(&(*begin));
        ++begin;
    }
}
//===================================================================================================
template <typename ElemT> class matrix_buffer_t {
    protected:
    ElemT* buf_  = nullptr;
    size_t rows_ = 0;
    size_t cols_ = 0;
    size_t used_ = 0;
//===================================================================================================
    protected:
    matrix_buffer_t(size_t rows, size_t cols): 
    buf_{static_cast<ElemT*>(::operator new(sizeof(ElemT) * rows * cols))}, rows_{rows}, cols_(cols) {}
    
    matrix_buffer_t& operator=(const matrix_buffer_t& other_buf) = delete;

    matrix_buffer_t(const matrix_buffer_t& other_buf)            = delete;

    matrix_buffer_t(matrix_buffer_t&& other_buf) noexcept: buf_ {other_buf.buf_},
                                                           rows_{other_buf.rows_},
                                                           cols_{other_buf.cols_},
                                                           used_{other_buf.used_} {
        other_buf.buf_ = nullptr;                                                           
        other_buf.rows_ = 0;
        other_buf.cols_ = 0;
        other_buf.used_ = 0;
    }

    matrix_buffer_t& operator=(matrix_buffer_t&& other_buf) noexcept {
        if (this == &other_buf)
            return *this; 

        std::swap(buf_,  other_buf.buf_ );
        std::swap(rows_, other_buf.rows_);
        std::swap(cols_, other_buf.cols_);
        std::swap(used_, other_buf.used_);

        return *this;
    }

    ~matrix_buffer_t() {
        Destructor(buf_, buf_ + used_);
        ::operator delete(buf_);
    }
};

} // namespace Matrix


#endif // BUFFER_HPP