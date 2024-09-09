#ifndef BUFFER_HPP
#define BUFFER_HPP

#include <algorithm>

namespace Matrix {
    template <typename T> void Constructor(T* p, const T& rhs) {
        new (p) T(rhs);
    }

    template <class T> void Destructor(T* p) {
        p->~T();
    }

    template <typename Iterator> void Destructor(Iterator begin, Iterator end) {
        while (begin != end) {
            Destructor(begin);
            ++begin;
        }
    }
//-------------------------------------------------------------------------------------
    template <typename T> class Buffer {
        protected:
        long long rows_ = 0;
        long long cols_ = 0;
        T*        buf_  = nullptr;
//-------------------------------------------------------------------------------------
        protected:
        #if 1
        Buffer(long long rows, long long cols): rows_(rows), cols_(cols) {
            if ((rows_ == 0) || (cols_ == 0))
                buf_ = nullptr;
            else 
                buf_ = static_cast<T*>(::operator new(sizeof(T) 
                                            * static_cast<unsigned long long>(rows) 
                                            * static_cast<unsigned long long>(cols)));
        }
        #endif
    
        Buffer(const Buffer<T>& other_buf) = delete;

        Buffer<T>& operator=(const Buffer<T>& other_buf) = delete;

        Buffer(Buffer<T>&& other) noexcept: rows_{other.rows_}
                                          , cols_{other.cols_}
                                          , buf_ {other.buf_ } {}

        Buffer<T>& operator=(Buffer<T>&& other_buf) noexcept {
            if (this == &other_buf)
                return *this; 

            std::swap(this->rows_, other_buf.rows_);                
            std::swap(this->cols_, other_buf.cols_);
            std::swap(this->buf_,  other_buf.buf_ );            
            return *this;
        }

        ~Buffer() {
            Destructor(buf_, buf_ + rows_ * cols_);
            ::operator delete(buf_);
        }
    };

} // namespace Matrix


#endif // BUFFER_HPP