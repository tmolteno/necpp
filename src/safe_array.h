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
#include <Eigen/Dense>

#include "nec_exception.h"

/* BoundsViol is always defined so that the throw target is consistent
   across translation units, regardless of NEC_ERROR_CHECK. */
class BoundsViol : public nec_exception {
public:
  BoundsViol(const char* message, int64_t index, int64_t bound)
    : nec_exception(message)
  {
    m_message << "array index: " << index << " exceeds " << bound << std::endl;
  }
};

/*! \brief A Safe Array class backed by Eigen for SIMD-accelerated operations.
 *
 * Bounds checking is done if the macro NEC_ERROR_CHECK is defined at compile time.
 * Storage is column-major (Eigen default) to match FORTRAN/NEC conventions.
 * Non-owning views (from segment()) use raw pointers to avoid Eigen expression
 * template lifetime issues.
 */
template<typename T>
class safe_array {
public:
  using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  safe_array()
    : _len(0), _rows(0), _cols(0), _capacity(0), _resize_chunk(2),
      _view_ptr(nullptr), _own_data(true)
  { }

  safe_array(int64_t in_size)
    : _len(0), _rows(0), _cols(0), _capacity(0), _resize_chunk(2),
      _view_ptr(nullptr), _own_data(true)
  {
    resize(in_size);
  }

  safe_array(const safe_array<T>& in_array)
    : _len(0), _rows(0), _cols(0), _capacity(0), _resize_chunk(2),
      _view_ptr(nullptr), _own_data(true)
  {
    copy(in_array);
  }

  ~safe_array() { /* Eigen::Matrix destructor handles cleanup */ }

  int64_t size() const    { return _len; }
  int64_t rows() const    { return _rows; }
  int64_t cols() const    { return _cols; }
  int64_t capacity() const { return _own_data ? _capacity : _len; }

  void resize(int64_t n_rows, int64_t n_cols) {
    _rows = n_rows;
    _cols = n_cols;
    resize(_rows * _cols);
  }

  void copy(const safe_array<T>& in_array) {
    if (in_array._rows == 0)
      resize(in_array._len);
    else
      resize(in_array._rows, in_array._cols);
    _storage.head(_len) = in_array._eigen_view();
  }

  void resize(int64_t new_length) {
#ifdef NEC_ERROR_CHECK
    if (!_own_data)
      throw new nec_exception("attempt to resize data we do not own");
#endif
    if (new_length > _capacity) {
      _capacity = new_length + new_length / 2;  // 1.5x growth
      try {
        Vector new_storage(_capacity);
        if (_len > 0)
          new_storage.head(_len) = _storage.head(_len);
        _storage.swap(new_storage);
      } catch (std::bad_alloc& ba) {
        throw new nec_exception("Error: Out of Memory ");
      }
    }
    _len = new_length;
  }

  T maxCoeff() const {
    if (0 == _len)
      throw new nec_exception("No elements in maxCoeff");
    return _eigen_view().maxCoeff();
  }

  T minCoeff() const {
    if (0 == _len)
      throw new nec_exception("No elements in minCoeff");
    return _eigen_view().minCoeff();
  }

  T sum(int64_t start_index, int64_t stop_index) {
    return _eigen_view().segment(start_index, stop_index - start_index).sum();
  }

  T sum() { return _eigen_view().sum(); }

  void fill(int64_t start, int64_t N, const T& x) {
    _eigen_view().segment(start, N).setConstant(x);
  }

  void setConstant(const T& x) { fill(0, _len, x); }

  void set_col_major(int64_t col_dim, int64_t col, int64_t row, const T& val) {
    (*this)[check(row * col_dim + col)] = val;
  }

  T& get_col_major(int64_t col_dim, int64_t col, int64_t row) {
    return (*this)[check(row * col_dim + col)];
  }

  T& getItem(int64_t row, int64_t col) {
    return (*this)[check(row, col)];
  }

  T& operator()(int64_t row, int64_t col) {
    return getItem(row, col);
  }

  const T& operator()(int64_t row, int64_t col) const {
    return (*this)[check(row, col)];
  }

  T& getItem(int64_t i) {
    return (*this)[check(i)];
  }

  const T& operator[](int64_t i) const {
    return _own_data ? _storage(check(i)) : _view_ptr[check(i)];
  }

  T& operator[](int64_t i) {
    return _own_data ? _storage(check(i)) : _view_ptr[check(i)];
  }

  safe_array<T> segment(int64_t start_index, int64_t end_index) {
    if (int64_t(-1) == end_index)
      throw "foo";
    int64_t n = end_index - start_index + 1;
    safe_array<T> result;
    result._view_ptr = data() + start_index;
    result._len = n;
    result._rows = 0;
    result._cols = 0;
    result._own_data = false;
    return result;
  }

  safe_array<T> eigen_segment(int64_t start_index, int64_t n) {
    return segment(start_index, start_index + n - 1);
  }

  T* data() {
    return _own_data ? _storage.data() : _view_ptr;
  }

  const T* data() const {
    return _own_data ? _storage.data() : _view_ptr;
  }

  safe_array<T>& operator=(const safe_array<T>& in_array) {
    copy(in_array);
    return *this;
  }

protected:
  int64_t _len;
  int64_t _rows;
  int64_t _cols;
  int64_t _capacity;
  int64_t _resize_chunk;

  Vector _storage;            // owning storage (unused when !_own_data)
  T*    _view_ptr;            // non-owning view pointer
  bool  _own_data;

  safe_array(const safe_array<T>& in_array, int64_t start_index,
             int64_t end_index, bool in_copy_data)
  {
    _resize_chunk = in_array._resize_chunk;
    _len = (end_index - start_index) + 1;
    _rows = 0;
    _cols = 0;

    if (in_copy_data) {
      _capacity = _len;
      _storage.resize(_len);
      _storage = in_array._eigen_view().segment(start_index, _len);
      _own_data = true;
      _view_ptr = nullptr;
    } else {
      _view_ptr = in_array.data() + start_index;
      _capacity = 0;
      _own_data = false;
    }
  }

  inline int64_t check(int64_t i) const {
#ifdef NEC_ERROR_CHECK
    if (i < 0 || i >= _len)
      throw new BoundsViol("safe_array: ", i, _len);
#endif
    return i;
  }

  inline int64_t check(int64_t row, int64_t col) const {
#ifdef NEC_ERROR_CHECK
    if (row < 0 || row >= _rows)
      throw new BoundsViol("safe_array: ", row, _rows);
    if (col < 0 || col >= _cols)
      throw new BoundsViol("safe_array: ", col, _cols);
#endif
    return check(int64_t(col) * _rows + row);
  }

private:
  Eigen::Map<Vector> _eigen_view() {
    return Eigen::Map<Vector>(data(), _len);
  }

  const Eigen::Map<const Vector> _eigen_view() const {
    return Eigen::Map<const Vector>(data(), _len);
  }
};


template<typename T>
class safe_matrix : public safe_array<T> {
public:
  safe_matrix(int64_t in_rows, int64_t in_cols) : safe_array<T>(in_rows * in_cols) {
    this->resize(in_rows, in_cols);
  }
  safe_matrix() : safe_array<T>() { }

private:
  using safe_array<T>::operator[];
};

#endif /* __safe_array__ */
