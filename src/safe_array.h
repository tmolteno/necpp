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
  
  You should have received a copy of the GNU General Public License with 
  this program; if not, write to the Free Software
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
    : _len(0), _rows(0), _cols(0), _resize_chunk(2), _data(NULL), _data_size(0), _own_data(true)
  { }

  safe_array(int64_t in_size)
    : _len(0), _rows(0), _cols(0), _resize_chunk(2), _data(NULL), _data_size(0), _own_data(true)
  {
    resize(in_size);
  }
  
  safe_array(const safe_array<T>& in_array)
    : _len(0), _rows(0), _cols(0), _resize_chunk(2), _data(NULL), _data_size(0), _own_data(true)
  {
    copy(in_array);
  }
  
  ~safe_array()  {
    if (_own_data) {
      delete[] _data;
      _data = NULL;
    }
  }
  
  
  int64_t size() const {
    return _len;
  }

  int64_t rows() const {
    return _rows;
  }

  int64_t cols() const {
    return _cols;
  }

  int64_t capacity() const {
    return _data_size;
  }

  void resize(int64_t n_rows, int64_t n_cols) {
    _rows = n_rows;
    _cols = n_cols;
    
    resize(_rows * _cols);
  }

  // Copy the contents of in_array to our array
  // resizing as appropriate.
  void copy(const safe_array<T>& in_array) {
    if (in_array._rows == 0)
      resize(in_array._len);
    else
      resize(in_array._rows,in_array._cols);
      
    for (int64_t i=0; i<_len; i++)
      _data[i] = in_array[i];
  }

  
  void resize(int64_t new_length) {
#ifdef NEC_ERROR_CHECK
    if (! _own_data)
      throw new nec_exception("attempt to resize data we do not own");
#endif
    if (new_length > _data_size) {
      // We allocate _resize_chunk more bytes than we need to avoid
      // resizing too often. 
      _data_size = new_length + _resize_chunk;
      try {
        T* new_data_ = new T[_data_size];
        if (0 != _len)
          std::memcpy(new_data_, _data, _len * sizeof(T));

        delete[] _data;
        _data = new_data_;
      } catch (std::bad_alloc& ba) {
        throw new nec_exception("Error: Out of Memory ");
      }
    }
    _len = new_length;
  }
  
  /*!\brief return the largest element of the array */
  T maxCoeff() const {
    if (0 == _len)
      throw new nec_exception("No elements in maxCoeff");

    T ret = _data[check(0)];
    
    for (int64_t i = 1; i < _len; i++ ) {
      if ( _data[check(i)] > ret)
        ret = _data[check(i)];
    }
    return ret;
  }

  /*!\brief return the largest element of the array */
  T minCoeff() const {
    if (int64_t(0) == _len)
      throw new nec_exception("No elements in minCoeff");

    T ret = _data[check(0)];
    
    for (int64_t i = 1; i < _len; i++ ) {
      if ( _data[check(i)] < ret)
        ret = _data[check(i)];
    }
    return ret;
  }

  /*!\brief return the sum of the specified elements in the array */
  T sum(int64_t start_index, int64_t stop_index)  {
    T ret = _data[check(start_index)];
    
    for (int64_t i = start_index+1; i < stop_index; i++ )  {
      ret += _data[check(i)];
    }
    return ret;
  }

  /*!\brief return the sum of all elements in the array */
  T sum() {
    return sum(0,_len);
  }

  // fill specified elements of the array with x
  void fill(int64_t start, int64_t N, const T& x) {
    int64_t stop = start + N;
    for (int64_t i = start; i < stop; i++ ) {
      _data[check(i)] = x;
    }
  }
  
  // fill all elements of the array with x
  void setConstant(const T& x)  {
    fill(0,_len,x);
  }

  /*! \brief Set an element assuming that the data is stored in column major form.
  */
  void set_col_major(int64_t col_dim, int64_t col, int64_t row, const T& val) {
    _data[check(row*col_dim + col)] = val;
  }

  /*! \brief Get an element assuming that the data is stored in column major form.
  */
  T& get_col_major(int64_t col_dim, int64_t col, int64_t row) {
    return _data[check(row*col_dim + col)];
  }
  
  
  /**
   * \remark This is an accessor function that is useful for wrapping.
   * */
  T& getItem(int64_t row, int64_t col) {
    return _data[check(row,col)];
  }

  
  /**
   * \remark We use round brackets for indexing into 2D arrays.
   * */
  T& operator()(int64_t row, int64_t col)  {
    return getItem(row, col);
  }
  
  const T& operator()(int64_t row, int64_t col) const  {
    return _data[check(row,col)];
  }
  
  /* \brief An accessor method to help with wrapping the C++ objects
  * into other languages (like python, ruby)
  */
  T& getItem(int64_t i) {
    return _data[check(i)];
  }
  
  const T& operator[](int64_t i) const  {
    return _data[check(i)];
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
    return _data;
  }
                  
  safe_array<T>& operator=(const safe_array<T>& in_array)  {
    copy(in_array);
    return *this;
  }
          

protected:
  int64_t _len;
  int64_t _rows;
  int64_t _cols;
  int64_t _resize_chunk;
  
  T*  _data;
  int64_t _data_size;
  
  bool _own_data;
  
  /*!\brief Constructor only used to construct segment
  \param start_index The first element of in_array to copy.
  \param end_index The last element of in_array to copy.
  \param in_copy_data True if we create a copy of data from in_array. False - just reference the data in in_array.
  */
  safe_array(const safe_array<T>& in_array, int64_t start_index, int64_t end_index, bool in_copy_data)
  {
    _resize_chunk = in_array._resize_chunk;
    _len = (end_index - start_index)+1; // include the end_index element
    _rows = 0;
    _cols = 0;
    
    if (in_copy_data)  {
      _data = new T[_len];
      _data_size = _len;
      
      for (int64_t i=0; i<_len; i++)
        _data[check(i)] = in_array[start_index + i];

      _own_data = true;
    } else {
      _data = in_array.data() + start_index;
      _data_size = 0;
      _own_data = false;
    }
  }
  
  inline int64_t check(int64_t i) const  {
#ifdef NEC_ERROR_CHECK
    if (i < 0 || i >= _len)
      throw new BoundsViol("safe_array: ", i, _len);
#endif
    return i;
  }

  inline int64_t check(int64_t row, int64_t col) const  {
#ifdef NEC_ERROR_CHECK
    if (row < 0 || row >= _rows)
      throw new BoundsViol("safe_array: ", row, _rows);
    if (col < 0 || col >= _cols)
      throw new BoundsViol("safe_array: ", col, _cols);
#endif
    return check(int64_t(col)*_rows + row);
  }
}; 


template<typename T>
class safe_matrix : public safe_array<T> {
public:
  safe_matrix(int64_t in_rows, int64_t in_cols) : safe_array<T>(in_rows*in_cols) {
    this->resize(in_rows, in_cols);
  }
  safe_matrix() : safe_array<T>()
  { }


private:
  using safe_array<T>::operator[]; /* Stop folk from using square bracket indexing */

};

#endif /* __safe_array__ */
