/***************************************************************************
 *   Copyright (C) 2004-2008 by Tim Molteno                                *
 *   tim@physics.otago.ac.nz                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#pragma once

#include <string>
#include <sstream>

//#include "common.h"

class nec_exception
{
public:
	nec_exception()
	{
	}

	nec_exception(const char* message)
	{
		m_message << message;
	}

	nec_exception(const char* message, int code)
	{
		m_message << message << code;
	}

	nec_exception(const nec_exception& other)
	{
		m_message << other.m_message.str();
	}

	template <class T> void append(const T& message)
	{
		m_message << message;
	}

	std::string get_message() const
	{
		return m_message.str();
	}

	static std::string string_printf(const char* fmt, ...);
	
protected:
	std::stringstream m_message;
};

#ifdef _MSC_VER
/*
	Visual C++ does not allow macros with variable argument lists. Therefore error messages
	will be meaningless when this is compiled using VC++, however at least it will compile!
*/
inline void nec_stop(const char* __fmt, ...)
{
	 nec_exception __nex("Undefined Error");
	// __nex->os_printf(__fmt, __VA_ARGS__);
	 throw __nex;
}
#else
#define nec_stop(__fmt, ...)\
{    nec_exception __nex;\
	std::string _mess = nec_exception::string_printf(__fmt, __VA_ARGS__); \
	__nex.append(_mess.c_str()); \
	throw __nex; \
}
#endif
