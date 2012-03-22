#ifndef __common__
#define __common__
/*
	Various Definitions for nec2++
	
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

/* Version information */
#define	nec_version_major "1"
#define	nec_version_minor "2"
#define	nec_version_build "3"

#ifndef nec_build_date
	#define nec_build_date "2005-07-29"
#endif

#ifndef build_version
	#define	nec_version nec_version_major "." nec_version_minor "." nec_version_build " [" nec_build_date "]"
#else
	#define nec_version build_version " [" nec_build_date "]"
#endif

/*
	These are some common constants that should be moved into more appropriate locations
*/


#define	ETA	376.73
#define	CVEL	299.8
#define	RETA	2.654420938E-3
#define	TOSP	1.128379167
#define ACCS	1.E-12
#define	SP	1.772453851
#define	FPI	12.56637062
#define	CONST2	4.771341188


#define	SMIN	1.e-3

#endif /* __common__ */
