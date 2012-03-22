/*
	Copyright (C) 2004-2005  Timothy C.A. Molteno
	
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
#include "nec_output.h"
#include "nec_exception.h"
#include <string>
#include <sstream>


/* ---------------------------------------------------------------------*/

nec_output_file::nec_output_file()
{
	set_file(NULL);
	set_error_mode(false);
}

void nec_output_file::set_file(FILE* in_fp)
{
	m_output_fp = in_fp;
	set_indent(0);
}

void nec_output_file::set_error_mode(bool f)
{
	m_error_mode = f;
}

/* private */
void nec_output_file::do_output(const char* str)
{
	if (NULL == m_output_fp)
		return;
	
	fprintf(m_output_fp, str);
	if (m_error_mode)
		fprintf(stderr,str);
}

void nec_output_file::endl(int n_lines)
{
	for (int i=0; i < n_lines; i++)
		do_output("\n");
		
	m_require_indent = true;
}

void nec_output_file::end_section()
{
	endl(3);
}

void nec_output_file::set_indent(int n)
{
	m_indent = n;
	m_require_indent = true;
	indent();
}

void nec_output_file::indent()
{
	if (m_require_indent)
	{
		for (int i=0; i< m_indent; i++)
			do_output(" ");
		
		m_require_indent = false;
	}
}

void nec_output_file::line(const char* in_str)
{
	string(in_str,true);
}

void nec_output_file::string(const char* in_str, bool require_endl)
{
	indent();
	do_output(in_str);
	if (require_endl)
		endl();
}

void nec_output_file::real(nec_float in_nec_float)
{
	real_out(11,4,in_nec_float,true);
}

void nec_output_file::integer(long in_integer)
{
	if (NULL == m_output_fp)
		return;
	
	fprintf(m_output_fp,"%ld",in_integer);
	if (m_error_mode)
		fprintf(stderr,"%ld",in_integer);
}

void nec_output_file::real_out(int w, int p, nec_float f, bool sci)
{
	if (NULL == m_output_fp)
		return;
	
	std::stringstream ss;
	ss << "%" << w << "." << p;
	
	if (sci)
		ss << "E";
	else
		ss << "f";
	
	std::string s = ss.str();
	const char* fmt = s.c_str();
	
	fprintf(m_output_fp,fmt,f);
	if (m_error_mode)
		fprintf(stderr,fmt,f);
}


#include "stdarg.h"
#include "safe_array.h"

void nec_output_file::nec_printf(const char* fmt, ...)
{
	if (NULL == m_output_fp)
		return;
	
	{	
	va_list ap; 		/* special type for variable    */
	
	safe_array<char> format(2048);	/* argument lists               */
	int count = 0;
	int i, j;		/* Need all these to store      */
	char c;			/* values below in switch       */
	double d;
	unsigned u;
	char *s;
	void *v;
	
	va_start(ap, fmt); 	/* must be called before work   */
	while (*fmt)
	{
		for (j = 0; fmt[j] && fmt[j] != '%'; j++)
			format[j] = fmt[j];                    /* not a format string          */
		if (j)
		{
			format[j] = '\0';
			count += fprintf(m_output_fp, format.data());    /* log it verbatim              */
			fmt += j;
		} 
		else
		{
			for (j = 0; !isalpha(fmt[j]); j++)
			{   /* find end of format specifier */
				format[j] = fmt[j];
				if (j && fmt[j] == '%')              /* special case printing '%'    */
					break;
			}
			format[j] = fmt[j];                    /* finish writing specifier     */
			format[j + 1] = '\0';                  /* don't forget NULL terminator */
			fmt += j + 1;
			
			switch (format[j])
			{                   /* cases for all specifiers     */
				case 'd':
				case 'i':                              /* many use identical actions   */
					i = va_arg(ap, int);                 /* process the argument         */
					count += fprintf(m_output_fp, format.data(), i); /* and log it                 */
					break;
				case 'o':
				case 'x':
				case 'X':
				case 'u':
					u = va_arg(ap, unsigned);
					count += fprintf(m_output_fp, format.data(), u);
					break;
				case 'c':
					c = (char) va_arg(ap, int);          /* must cast!                   */
					count += fprintf(m_output_fp, format.data(), c);
					break;
				case 's':
					s = va_arg(ap, char *);
					count += fprintf(m_output_fp, format.data(), s);
					break;
				case 'f':
				case 'e':
				case 'E':
				case 'g':
				case 'G':
					d = va_arg(ap, double);
					count += fprintf(m_output_fp, format.data(), d);
					break;
				case 'p':
					v = va_arg(ap, void *);
					count += fprintf(m_output_fp, format.data(), v);
					break;
				case 'n':
					count += fprintf(m_output_fp, "%d", count);
					break;
				case '%':
					count += fprintf(m_output_fp, "%%");
					break;
				default:
					throw new nec_exception("Invalid format specifier in nec_printf()");
			}
		}
	}
	
	va_end(ap);	/* clean up */
	}
}

