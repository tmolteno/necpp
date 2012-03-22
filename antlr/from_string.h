/*
	Utility file for parsing text in C++
	
	Copyright (C) 2008  Timothy C.A. Molteno (tim@physics.otago.ac.nz)
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __from_string__
#define __from_string__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

template <class T> T from_string(const std::string& s)
{
        T t;
        std::istringstream iss(s);
        iss >> std::dec >> t;
        if (iss.fail())
        {
                std::cout << "Failed to parse " << s << std::endl;
                throw -1;
        }

        return t;
}


#endif /* __from_string__ */
