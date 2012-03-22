#ifndef __nec_output__
#define __nec_output__

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

#include <stdio.h>
#include "math_util.h"
#include "nec_debug.h"

class nec_output_flags
{
public:
	nec_output_flags()
	{
		do_gain = false;
		do_nec_output = true;
	}

	void set_gain_only(bool in_gain_only)
	{
		do_gain = in_gain_only;
		do_nec_output = !in_gain_only;
	}
	
	bool get_nec_flag()
	{
		return do_nec_output;
	}
	
	bool get_gain_flag()
	{
		return do_gain;
	}
	
private:	
	bool do_gain;
	bool do_nec_output;
};

class nec_output_file
{
public:
	nec_output_file();
	
	void set_file(FILE* in_fp);

	void endl(int n_lines = 1);

	void end_section();

	void set_indent(int n = 0);

	void line(const char* in_str);
	void string(const char* in_str, bool require_endl = false);
	void real(nec_float in_nec_float);
	void real_out(int w, int p, nec_float f, bool sci = true);
	void integer(long in_integer);
	
	void set_error_mode(bool f);
	
	void nec_printf(const char* fmt, ...);
					   	
	FILE* get_fp()
	{
		return m_output_fp;
	}
	
private:
	void indent();
	void do_output(const char* str);

	FILE* m_output_fp;
	bool m_require_indent;
	int m_indent;
	
	bool m_error_mode;
	FILE* m_error_fp;
};


/**
	A little class for setting up error mode, and then
	automatically going back to normal output mode.
	Usage
	
	{
		nec_error_mode em(s_output);
		s_output.line("Darn...");
	}
	
	s_output.line("Normal output here");
*/
class nec_error_mode
{
public:
	nec_error_mode(nec_output_file& of)
		: m_of(of)
	{
		of.set_error_mode(true);
	}
	
	~nec_error_mode()
	{
		m_of.set_error_mode(false);
	}
private:
	nec_output_file& m_of;	
};

#endif /* __nec_output__ */
