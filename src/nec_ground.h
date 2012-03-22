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
#ifndef __nec_ground__
#define __nec_ground__

#include "common.h"
#include "math_util.h"
#include "nec_debug.h"

/*
	\section{NEC Ground Specification}
	
	\subsection{the GN Card}
	
	The GN card specifies the relative dielectric constant and conductivity of ground in the vicinity of the antenna. 
	In addition, a second set of ground parameters for a second medium can be specified,
	or a radial wire ground screen can be modeled using a reflection coefficient approximation.
	
	The parameters of the second medium can also be specified on another data card whose mnemonic is GD.
	With the GD card, the parameters of the second medium can be varied and only the radiated fields
	need to be recalculated.
	Furthermore, if a radial wire ground screen has been specified on the GN card,
	the GD card is the only way to include a second medium.
	See Section~\ref{sec_card_gd} for details.
	
	\subsection{"gd" card: ground representation}
	\label{sec_card_gd}
	
	To specify the ground parameters of a second medium which is not in the immediate vicinity of the antenna.
	This card may only be used if a GN card has also been used.
	It does not affect the fields of surface patches.
	
	<ground type=finite>
		<dielectric_constant = 80/>
		<conductivity = 4/>
		<radial-wire>
			<count=12>
			<length=15 unit=cm>
			<
		</radial-wire>
	</ground>
	
*/
class nec_ground
{
public:

	nec_ground()
	{
		default_values();
	}
	
	nec_ground(const nec_ground& in_ground)
	{
		iperf = in_ground.iperf;
		ksymp = in_ground.ksymp;
		epsr = in_ground.epsr;	// relative dielectric constant
		sig = in_ground.sig;	// Conductivity
		
		// radial wire ground
		radial_wire_count = in_ground.radial_wire_count;
		radial_wire_length = in_ground.radial_wire_length;
		radial_wire_radius = in_ground.radial_wire_radius;
		
		// second medium parameters
		cliff_edge_distance = in_ground.cliff_edge_distance;
		cliff_height = in_ground.cliff_height;
		epsr2 = in_ground.epsr2;	// Relative dielectric constant
		sig2 = in_ground.sig2;		// Conductivity in mhos/meter
	}
	
	void default_values()
	{
		zrati=cplx_10();
		ksymp=1;
		radial_wire_count = 0;
		radial_wire_length = 0.0;
		radial_wire_radius = 0.0;
		iperf=0;
		/*22/09/2005 : added initialization for sig and epsr*/
		epsr = 0;
		sig = 0;
	}
	/*
		Parse a GN card. The input parameters here are the fields of the
		GN card.	
	*/
	void parse_gn(int itmp1, int itmp2,
		nec_float tmp1, nec_float tmp2,
		nec_float tmp3, nec_float tmp4,
		nec_float tmp5, nec_float tmp6
		);
	
	/*
		Setup a cliff (two medium ground)
	*/
	void setup_cliff(nec_float in_eprs2,
		nec_float in_sig2,
		nec_float clt, nec_float cht);
	
	nec_complex get_zrati2(nec_float wavelength);
	
	nec_float get_cl(nec_float wavelength) // cliff edge in wavelengths, 
	{
		return cliff_edge_distance / wavelength;
	}
	
	nec_float get_ch(nec_float wavelength) // cliff Height in wavelengths.
	{
		return cliff_height / wavelength;
	}
	
	nec_complex get_zrati_sqr() const
	{
		return zrati*zrati;
	}
	
	// accessors for the ground type
	inline bool type_finite_reflection()	{	return (0 == iperf); }
	inline bool type_perfect()				{	return (1 == iperf); }
	inline bool type_sommerfeld_norton()	{	return (2 == iperf); }
	
	
	bool is_valid() const
	{
		if (iperf < 0) return false;
		if (iperf > 2) return false;
		
		if (ksymp < 1) return false;
		if (ksymp > 2) return false;
		
		return true;
	}
	
	/*! \brief Return true if a ground is present */
	bool present() const
	{
		ASSERT(is_valid());
		
		if (2 == ksymp)
			return true;
			
		return false;
	}
	
	nec_float epsr;	// relative dielectric constant
	nec_float sig;	// Conductivity
	
	// radial wire ground
	int radial_wire_count;
	nec_float radial_wire_length;
	nec_float radial_wire_radius;
	
	// second medium parameters
	nec_float cliff_edge_distance;
	nec_float cliff_height;
	nec_float epsr2;	// Relative dielectric constant
	nec_float sig2;		// Conductivity in mhos/meter

//	MOVE THESE TO THE GROUND
	
	nec_float scrwl;	// length of wires in radial-wire ground-screen in wavelengths
	nec_float scrwr;	// radius of wires in radial-wire ground-screen in wavelengths
	
	nec_complex m_t1;
	nec_float t2;	// constants for the radial-wire ground-screen impedance
	
	nec_complex zrati, frati;
	
	nec_float t1;	// constants for the radial-wire ground-screen impedance
	
	/** Ground Flag. ==1 no ground, =2 ground present */
	int ksymp;
private:
	/* iperf: Ground-type flag. The options are: 
		-1 - nullifies ground parameters previously used and sets free- space condition. The remainder of the card is left blank in this case. 
		O - finite ground, reflection-coefficient approximation. 
		1 - perfectly conducting ground. 
		2 - finite ground, Sommerfeld/Norton method.
	*/
	int iperf;

	
};

#endif /* __nec_ground__ */

