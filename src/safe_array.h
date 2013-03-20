/*
	Copyright (C) 2004-2013  Timothy C.A. Molteno
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __safe_array__
#define __safe_array__

#include <iostream>
#include <cstring>
#include <sstream>

#include "nec_exception.h"

#ifdef NEC_ERROR_CHECK


class BoundsViol : public nec_exception
{
public:
	BoundsViol(const char* message, long index, long bound)
		: nec_exception(message)
	{
		m_message << "array index: " << index << " exceeds " << bound << std::endl;
	}
};
#endif

/*! \brief 	A Safe Array class for nec2++ that performs bounds checking if the macro NEC_ERROR_CHECK is defined at compile time.
	
	This class also includes some utility functions for handling
	common vector operations.
	
	\todo Modify the resize() operator so  that we resize in units of _resize_chunk. This will
	increase efficiency.
*/
template<typename T>
class safe_array
{
public:
	safe_array()
		: len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
	{ }

	safe_array(long in_size)
		: len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
	{
		resize(in_size);
	}
	
/*	safe_array(long n_rows, long n_cols)
		: len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
	{
		resize(n_rows, n_cols);
	}*/
	
	safe_array(const safe_array<T>& in_array)
		: len_(0), rows_(0), cols_(0), resize_chunk_(2), data_(NULL), data_size_(0), own_data_(true)
	{
		copy(in_array);
	}
	
	~safe_array()
	{
		if (own_data_)
			delete[] data_;
	}
	
	
	long size() const
	{
		return len_;
	}

	void resize(long n_rows, long n_cols)
	{
		rows_ = n_rows;
		cols_ = n_cols;
		
		resize(rows_ * cols_);
	}

	// Copy the contents of in_array to our array
	// resizing as appropriate.
	void copy(const safe_array<T>& in_array)
	{
		if (in_array.rows_ == 0)
			resize(in_array.len_);
		else
			resize(in_array.rows_,in_array.cols_);
			
		for (long i=0; i<len_; i++)
			data_[i] = in_array[i];
	}

	
	void resize(long new_length)
	{
#ifdef NEC_ERROR_CHECK
		if (! own_data_)
			throw new nec_exception("attempt to resize data we do not own");
#endif			
		if (new_length > data_size_)
		{
			// We allocate resize_chunk_ more bytes than we need to avoid
			// resizing too often. 
			T* new_data_ = new T[new_length + resize_chunk_];
			
			data_size_ = new_length + resize_chunk_;
			
			if (0 != len_)
				std::memcpy(new_data_, data_, len_ * sizeof(T));
		
			delete[] data_;
			data_ = new_data_;
			len_ = new_length;
		}
		else
		{
			len_ = new_length;
		}
	}
	
	/*!\brief return the largest element of the array */
	T maxCoeff() const
	{
		if (0 == len_)
			throw new nec_exception("No elements in maxCoeff");

		T ret = data_[check(0)];
		
		for (long i = 1; i < len_; i++ )
		{
			if ( data_[check(i)] > ret)
				ret = data_[check(i)];
		}
		return ret;
	}

	/*!\brief return the largest element of the array */
	T minCoeff() const
	{
		if (0 == len_)
			throw new nec_exception("No elements in minCoeff");

		T ret = data_[check(0)];
		
		for (long i = 1; i < len_; i++ )
		{
			if ( data_[check(i)] < ret)
				ret = data_[check(i)];
		}
		return ret;
	}

	/*!\brief return the sum of the specified elements in the array */
	T sum(long start_index, long stop_index)
	{
		T ret = data_[check(start_index)];
		
		for (long i = start_index+1; i < stop_index; i++ )
		{
			ret += data_[check(i)];
		}
		return ret;
	}

	/*!\brief return the sum of all elements in the array */
	T sum()
	{
		return sum(0,len_);
	}

	// fill specified elements of the array with x
	void fill(long start, long N, const T& x)
	{
		long stop = start + N;
		for (long i = start; i < stop; i++ )
		{
			data_[check(i)] = x;
		}
	}
	
	// fill all elements of the array with x
	void setConstant(const T& x)
	{
		fill(0,len_,x);
	}

	/*! \brief Set an element assuming that the data is stored in column major form.
	*/
	void set_col_major(int col_dim, int col, int row, const T& val)
	{
		data_[check(row*col_dim + col)] = val;
	}

	/*! \brief Get an element assuming that the data is stored in column major form.
	*/
	T& get_col_major(int col_dim, int col, int row)
	{
		return data_[check(row*col_dim + col)];
	}
	
	T& operator()(long row, long col)
	{
		return data_[check(row,col)];
	}
	
	const T& operator[](long i) const
	{
		return data_[check(i)];
	}
	
	T& operator[](long i)
	{
		return data_[check(i)];
	}
	
	/*!\brief Return a representation of a subset of this array
		 \param start_index The index of the first element
		 \param end_index If -1, then finish at the end of the array
	*/
	safe_array<T> segment(long start_index, long end_index)
	{
		if (-1 == end_index)
			throw "foo";
		return eigen_segment(start_index, end_index);
/*		if (-1 == end_index)
			end_index = len_;
			
		return safe_array<T>(*this, start_index, end_index, false);*/
	}
	
	/*!\brief Return a representation of a subset of this array
		 \param start_index The index of the first element
		 \param end_index If -1, then finish at the end of the array
	*/
	safe_array<T> eigen_segment(long start_index, long n)
	{
		long end_index = start_index + n;
			
		return safe_array<T>(*this, start_index, end_index, false);
	}
	
	
	T* data() const
	{
		return data_;
	}
			
	safe_array<T>& operator=(const safe_array<T>& in_array)
	{
		copy(in_array);
		return *this;
	}
		

private:
	long len_;
	long rows_;
	long cols_;
	long resize_chunk_;
	
	T*  data_;
	long data_size_;
	
	bool own_data_;
	
	/*!\brief Constructor only used to construct segment
	\param in_copy_data True if we create a copy of data from in_array. False - just reference the data in in_array.
	*/
	safe_array(safe_array<T>& in_array, long start_index, long end_index, bool in_copy_data)
	{
		resize_chunk_ = in_array.resize_chunk_;
		len_ = (end_index - start_index);
		rows_ = 0;
		cols_ = 0;
		
		if (in_copy_data)
		{
			data_ = new T[len_];
			data_size_ = len_;
			
			for (long i=0; i<len_; i++)
				data_[check(i)] = in_array[start_index + i];
		
			own_data_ = true;
		}
		else
		{
			data_ = in_array.data() + start_index;
			data_size_ = 0;
			own_data_ = false;
		}
	}
	
	inline long check(long i) const
	{
#ifdef NEC_ERROR_CHECK
		if (i < 0 || i >= len_)
			throw new BoundsViol("safe_array: ", i, len_);
#endif
		return i;
	}

	inline long check(long row, long col) const
	{
#ifdef NEC_ERROR_CHECK
		if (row < 0 || row >= rows_)
			throw new BoundsViol("safe_array: ", row, rows_);
		if (col < 0 || col >= cols_)
			throw new BoundsViol("safe_array: ", col, cols_);
#endif
		return check(col*rows_ + row);
	}
}; 


#endif /* __safe_array__ */
