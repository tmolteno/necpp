/***************************************************************************
 *   Copyright (C) 2004 by Tim Molteno                                     *
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


#include "safe_array.h"
#include "nec_exception.h"
#include "stdarg.h"

using namespace std;

string nec_exception::string_printf(const char* fmt, ...)
{
	stringstream _sstream;
	va_list ap; 		/* special type for variable    */
	
	safe_array<char> format(2048);	/* argument lists               */
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
			_sstream << format.data();    /* log it verbatim              */
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
			
			switch (format[j]) {                   /* cases for all specifiers     */
			case 'd':
			case 'i':                              /* many use identical actions   */
				i = va_arg(ap, int);                 /* process the argument         */
				_sstream <<  i;							/* and log it                 */
				break;
			case 'o':
			case 'x':
			case 'X':
			case 'u':
				u = va_arg(ap, unsigned);
				_sstream << u;
				break;
			case 'c':
				c = (char) va_arg(ap, int);          /* must cast!                   */
				_sstream << c;
				break;
			case 's':
				s = va_arg(ap, char *);
				_sstream << s;
				break;
			case 'f':
			case 'e':
			case 'E':
			case 'g':
			case 'G':
				d = va_arg(ap, double);
				_sstream << d;
				break;
			case 'p':
				v = va_arg(ap, void *);
				_sstream << v;
				break;
			case '%':
				_sstream << "%%";
				break;
			default:
				throw new nec_exception("Invalid format specifier in os_printf()");
			}
		}
	}
	
	va_end(ap);	/* clean up */
	return _sstream.str();
}
