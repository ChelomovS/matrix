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
    int rows_   = 0;
    int cols_   = 0;
    ElemT* buf_ = nullptr;
//===================================================================================================
    protected:
    matrix_buffer_t(int rows, int cols): rows_{rows}, 
                                         cols_{cols} {
        buf_ = static_cast<ElemT*>(::operator new(sizeof(ElemT) * 
                                                  static_cast<unsigned int>(rows_) *
                                                  static_cast<unsigned int>(cols_)));
    }
    
    matrix_buffer_t& operator=(const matrix_buffer_t& other_buf) = delete;

    matrix_buffer_t(const matrix_buffer_t& other_buf)            = delete;

    matrix_buffer_t(matrix_buffer_t&& other_buf) noexcept: rows_{other_buf.rows_}, 
                                                           cols_{other_buf.cols_},
                                                           buf_ {other_buf.buf_} {
        other_buf.buf_ = nullptr;
    }

    matrix_buffer_t& operator=(matrix_buffer_t&& other_buf) noexcept {
        if (this == &other_buf)
            return *this; 

        std::swap(this->rows_, other_buf.rows_);                
        std::swap(this->cols_, other_buf.cols_);
        std::swap(this->buf_,  other_buf.buf_ );            
        return *this;
    }

    ~matrix_buffer_t() {
        ::operator delete(buf_);
    }
};

} // namespace Matrix


#endif // BUFFER_HPP