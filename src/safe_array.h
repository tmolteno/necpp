/*
  Copyright (C) 2004-2013,2015  Timothy C.A. Molteno
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  aint32_t with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __safe_array__
#define __safe_array__

#include <iostream>
#include <cstring>
#include <sstream>
#include <stdint.h>

#include "nec_exception.h"

#ifdef NEC_ERROR_CHECK


class BoundsViol : public nec_exception {
public:
  BoundsViol(const char* message, int64_t index, int64_t bound)
    : nec_exception(message)
  {
    m_message << "array index: " << index << " exceeds " << bound << std::endl;
  }
};
#endif

/*! \brief A Safe Array class for nec2++ that performs bounds checking 
 * 
 * Bounds checking is done if the macro NEC_ERROR_CHECK is defined at compile time.
 * 
 * This class also includes some utility functions for handling common vector operations.
 */
template<typename T>
class safe_array {
public:
  safe_array()
    : len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
  { }

  safe_array(int64_t in_size)
    : len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
  {
    resize(in_size);
  }
  
  safe_array(const safe_array<T>& in_array)
    : len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
  {
    copy(in_array);
  }
  
  ~safe_array()  {
    if (own_data_) {
      delete[] data_;
      data_ = NULL;
    }
  }
  
  
  int64_t size() const {
    return len_;
  }

  int32_t rows() const {
    return rows_;
  }

  int32_t cols() const {
    return cols_;
  }

  int64_t capacity() const {
    return data_size_;
  }

  void resize(int32_t n_rows, int32_t n_cols) {
    rows_ = n_rows;
    cols_ = n_cols;
    
    resize(rows_ * cols_);
  }

  // Copy the contents of in_array to our array
  // resizing as appropriate.
  void copy(const safe_array<T>& in_array) {
    if (in_array.rows_ == 0)
      resize(in_array.len_);
    else
      resize(in_array.rows_,in_array.cols_);
      
    for (int64_t i=0; i<len_; i++)
      data_[i] = in_array[i];
  }

  
  void resize(int64_t new_length) {
#ifdef NEC_ERROR_CHECK
    if (! own_data_)
      throw new nec_exception("attempt to resize data we do not own");
#endif
    if (new_length > data_size_) {
      // We allocate resize_chunk_ more bytes than we need to avoid
      // resizing too often. 
      data_size_ = new_length + resize_chunk_;
      try {
        T* new_data_ = new T[data_size_];
        if (0 != len_)
          std::memcpy(new_data_, data_, len_ * sizeof(T));

        delete[] data_;
        data_ = new_data_;
      } catch (std::bad_alloc& ba) {
        throw new nec_exception("Error: Out of Memory ");
      }
    }
    len_ = new_length;
  }
  
  /*!\brief return the largest element of the array */
  T maxCoeff() const {
    if (0 == len_)
      throw new nec_exception("No elements in maxCoeff");

    T ret = data_[check(0)];
    
    for (int64_t i = 1; i < len_; i++ ) {
      if ( data_[check(i)] > ret)
        ret = data_[check(i)];
    }
    return ret;
  }

  /*!\brief return the largest element of the array */
  T minCoeff() const {
    if (int64_t(0) == len_)
      throw new nec_exception("No elements in minCoeff");

    T ret = data_[check(0)];
    
    for (int64_t i = 1; i < len_; i++ ) {
      if ( data_[check(i)] < ret)
        ret = data_[check(i)];
    }
    return ret;
  }

  /*!\brief return the sum of the specified elements in the array */
  T sum(int64_t start_index, int64_t stop_index)  {
    T ret = data_[check(start_index)];
    
    for (int64_t i = start_index+1; i < stop_index; i++ )  {
      ret += data_[check(i)];
    }
    return ret;
  }

  /*!\brief return the sum of all elements in the array */
  T sum() {
    return sum(0,len_);
  }

  // fill specified elements of the array with x
  void fill(int64_t start, int64_t N, const T& x) {
    int64_t stop = start + N;
    for (int64_t i = start; i < stop; i++ ) {
      data_[check(i)] = x;
    }
  }
  
  // fill all elements of the array with x
  void setConstant(const T& x)  {
    fill(0,len_,x);
  }

  /*! \brief Set an element assuming that the data is stored in column major form.
  */
  void set_col_major(int32_t col_dim, int32_t col, int32_t row, const T& val) {
    data_[check(row*col_dim + col)] = val;
  }

  /*! \brief Get an element assuming that the data is stored in column major form.
  */
  T& get_col_major(int32_t col_dim, int32_t col, int32_t row) {
    return data_[check(row*col_dim + col)];
  }
  
  
  /**
   * \remark This is an accessor function that is useful for wrapping.
   * */
  T& getItem(int32_t row, int32_t col) {
    return data_[check(row,col)];
  }

  
  /**
   * \remark We use round brackets for indexing into 2D arrays.
   * */
  T& operator()(int32_t row, int32_t col)  {
    return getItem(row, col);
  }
  
  const T& operator()(int32_t row, int32_t col) const  {
    return data_[check(row,col)];
  }
  
  /* \brief An accessor method to help with wrapping the C++ objects
  * into other languages (like python, ruby)
  */
  T& getItem(int64_t i) {
    return data_[check(i)];
  }
  
  const T& operator[](int64_t i) const  {
    return data_[check(i)];
  }
  
  T& operator[](int64_t i)  {
    return getItem(i);
  }
  
  /*!\brief Return a representation of a subset of this array
      \param start_index The index of the first element
      \param end_index If -1, then finish at the end of the array
  */
  safe_array<T> segment(int64_t start_index, int64_t end_index)  {
    if (int64_t(-1) == end_index)
      throw "foo";
    return safe_array<T>(*this, start_index, end_index, false);
  }
  
  /*!\brief Return a representation of a subset of this array
      \param start_index The index of the first element
      \param n Number of elements in the segment
  */
  safe_array<T> eigen_segment(int64_t start_index, int64_t n)  {
    int64_t end_index = start_index + n;
    return safe_array<T>(*this, start_index, end_index, false);
  }
  
  
  T* data() const  {
    return data_;
  }
                  
  safe_array<T>& operator=(const safe_array<T>& in_array)  {
    copy(in_array);
    return *this;
  }
          

protected:
  int64_t len_;
  int32_t rows_;
  int32_t cols_;
  int64_t resize_chunk_;
  
  T*  data_;
  int64_t data_size_;
  
  bool own_data_;
  
  /*!\brief Constructor only used to construct segment
  \param start_index The first element of in_array to copy.
  \param end_index The last element of in_array to copy.
  \param in_copy_data True if we create a copy of data from in_array. False - just reference the data in in_array.
  */
  safe_array(const safe_array<T>& in_array, int64_t start_index, int64_t end_index, bool in_copy_data)
  {
    resize_chunk_ = in_array.resize_chunk_;
    len_ = (end_index - start_index)+1; // include the end_index element
    rows_ = 0;
    cols_ = 0;
    
    if (in_copy_data)  {
      data_ = new T[len_];
      data_size_ = len_;
      
      for (int64_t i=0; i<len_; i++)
        data_[check(i)] = in_array[start_index + i];

      own_data_ = true;
    } else {
      data_ = in_array.data() + start_index;
      data_size_ = 0;
      own_data_ = false;
    }
  }
  
  inline int64_t check(int64_t i) const  {
#ifdef NEC_ERROR_CHECK
    if (i < 0 || i >= len_)
      throw new BoundsViol("safe_array: ", i, len_);
#endif
    return i;
  }

  inline int64_t check(int32_t row, int32_t col) const  {
#ifdef NEC_ERROR_CHECK
    if (row < 0 || row >= rows_)
      throw new BoundsViol("safe_array: ", row, rows_);
    if (col < 0 || col >= cols_)
      throw new BoundsViol("safe_array: ", col, cols_);
#endif
    return check(int64_t(col)*rows_ + row);
  }
}; 


template<typename T>
class safe_matrix : public safe_array<T> {
public:
  safe_matrix(int32_t _rows, int32_t _cols) : safe_array<T>(_rows*_cols) {
    this->resize(_rows, _cols);
  }
  safe_matrix() : safe_array<T>()
  { }


private:
  using safe_array<T>::operator[]; /* Stop folk from using square bracket indexing */

};

#endif /* __safe_array__ */
