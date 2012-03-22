/***************************************************************************
 *   Copyright (C) 2004-2005 by Tim Molteno                                *
 *   tim@molteno.net                                                       *
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
#ifndef __nec_debug__
#define __nec_debug__

#include <iostream>
#include <string>
#include <sstream>

#ifdef NEC_ERROR_CHECK
	#define DEBUG_TRACE(__x) {std::cout << __x << std::endl;}
	#define ASSERT(__x) \
		{ if (false == (__x))\
			{	std::stringstream __ss; \
				__ss << "assert in file " << __FILE__ << " at line " << __LINE__;\
				std::string __s = __ss.str(); \
				throw __s.c_str(); \
			} \
		}
	#define ASSERT_EQUAL(__x, __y) ASSERT(fabs((__x) - (__y)) < 1e-15)
#else
	#define DEBUG_TRACE(__x)
	#define ASSERT(__x)
	#define ASSERT_EQUAL(__x, __y)
#endif /* NEC_ERROR_CHECK */

#endif /* __nec_debug__ */
