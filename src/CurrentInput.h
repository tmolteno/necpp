/*
	Copyright (C) 2004  Timothy C.A. Molteno
	
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
#ifndef __Current_Input__
#define __Current_Input__

#include "math_util.h"
#include <vector>

#include "BaseInput.h"

class segment
{
public:
	int number;
	int tag;
	nec_float x,y,z;
	nec_float length;

	nec_complex current;

	// Read the segment data from a NEC-2 output file
	segment(istream& m_stream)
	{
	}
};

class CurrentInput : public BaseInput
{
public:	
	vector<segment> segments;

	long n_items;
	
	CurrentInput(std::string& filename)
		: BaseInput(filename)
	{
		n_items = 0;
		
		string searchString("CURRENTS AND LOCATION");
		while (m_stream.good())
		{
			string line = readline();
	
			if (line.find(searchString,0) != string::npos)
			{
				
				while (line.find("PHASE",0) == string::npos)
					line = readline();
	
				line = readline();
				
				while (line != "")
				{
					stringstream ss(line);
				
					segment s(ss);
					segments.push_back(s);
	
					line = readline();
					n_items++;
				}
				cout << "Currents and Location: " << n_items << " lines" << endl;
			}	
		}
	}

	bool equalto(const CurrentInput& ai)
	{
		if (difference(ai) > 1e-5)
			return false;

		return true;
	}

	nec_float difference(const CurrentInput& ai)
	{
		nec_float ret = 0.0;
		
		if (n_items != ai.n_items)
			return 1;
		
		for (long i=0; i < n_items;i++)
		{
			try
			{
				ret += segments[i].diff(ai.segments[i]);
			}
			catch (string message)
			{
				cout << "Diff at segment [" << i << "] : " << message << endl;
			}
		}
		return  ret;
	};
};

#endif /* __Current_Input__ */
