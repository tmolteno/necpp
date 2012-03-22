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

	This class is not fully functional yet. It will perform
	the network solution when done!
*/
#include "matrix_algebra.h"
#include "electromag.h"
#include "nec_ground.h"
#include "c_geometry.h"

#include <stdio.h>

#include "nec_results.h"
extern nec_results s_results;

#include "nec_output.h"
extern nec_output_file s_output;
extern nec_output_flags s_output_flags;
extern nec_float wavelength;
extern int nop;
extern complex_array symmetry_array;


/* subroutine netwk solves for structure currents for a given */
/* excitation including the effect of non-radiating networks if */
/* present. */
class c_network
{
public:
	c_network();
	
	void net_solve(
		complex_array& cm,
		nec_complex *cmb,
		nec_complex *cmc,
		nec_complex *cmd,
		int_array& ip,
		complex_array& einc );
	
private:
	
	FILE *output_fp;
	
	c_geometry geometry;
	nec_ground ground;
	/* common  /crnt/
	
	*/
	real_array air, aii; // coefficients of the constant terms in the current interpolation functions for the current vector 
	real_array bir, bii; // coefficients of the sine terms in the current interpolation functions
	real_array cir, cii; // coefficients of the cosine terms in the current interpolation functions
	complex_array current_vector;	// the current vector
	
	/* common  /vsorc/ */
	int_array ivqd, source_segment_array, iqds;
	int nvqd, voltage_source_count, nqds;
	complex_array vqd, vqds, source_voltage_array;
	
	/* common  /netcx/ */
	int masym, neq, npeq, neq2, network_count, ntsol, nprint;
	int_array iseg1, iseg2, ntyp;
	real_array x11r, x11i, x12r;
	real_array x12i, x22r, x22i;
	nec_float input_power, network_power_loss;
	nec_complex zped;
};

void c_network::net_solve( complex_array& cm, nec_complex *cmb,
    nec_complex *cmc, nec_complex *cmd, int_array& ip,
    complex_array& einc )
{	
	/* Network buffers */
	int_array ipnt, nteqa, ntsca;
	complex_array vsrc, rhs, cmn, rhnt, rhnx;
	
	bool jump1, jump2;
	
	int nteq=0, ntsc=0, nseg2, irow2=0;
	int neqz2, neqt, irow1=0, i, nseg1, isc1=0, isc2=0;
	nec_float asmx, asa, y11r, y11i, y12r, y12i, y22r, y22i;
	nec_complex ymit, vlt, cux;
	
	neqz2= neq2;
	if ( neqz2 == 0)
		neqz2=1;
	
	input_power = 0.0;
	network_power_loss = 0.0;
	neqt= neq+ neq2;
	
	int ndimn = (2*network_count + voltage_source_count);
	
	/* Allocate network buffers */
	if ( network_count > 0 )
	{
		rhs.resize( geometry.n_plus_3m ); // this should probably be ndimn!
	
		rhnt.resize( ndimn );
		rhnx.resize( ndimn);
		cmn.resize( ndimn * ndimn );
	
		ntsca.resize( ndimn );
		nteqa.resize( ndimn );
		ipnt.resize( ndimn );
	
		vsrc.resize( voltage_source_count );
	}
	
	if ( ntsol == 0)
	{
		/* compute relative matrix asymmetry */
		if ( masym != 0)
		{
			irow1=0;
			for( i = 0; i < network_count; i++ )
			{
				nseg1= iseg1[i];
				for( isc1 = 0; isc1 < 2; isc1++ )
				{
					if ( irow1 == 0)
					{
						ipnt[irow1]= nseg1;
						nseg1= iseg2[i];
						irow1++;
						continue;
					}
			
					int j = 0;
					for( j = 0; j < irow1; j++ )
						if ( nseg1 == ipnt[j])
							break;
			
					if ( j == irow1 )
					{
						ipnt[irow1]= nseg1;
						irow1++;
					}
			
					nseg1= iseg2[i];
				} /* for( isc1 = 0; isc1 < 2; isc1++ ) */
			} /* for( i = 0; i < network_count; i++ ) */
	
			ASSERT(voltage_source_count >= 0);
			for( i = 0; i < voltage_source_count; i++ )
			{
				nseg1= source_segment_array[i];
				if ( irow1 == 0)
				{
					ipnt[irow1]= nseg1;
					irow1++;
					continue;
				}
			
				int j = 0;
				for( j = 0; j < irow1; j++ )
					if ( nseg1 == ipnt[j])
						break;
			
				if ( j == irow1 )
				{
					ipnt[irow1]= nseg1;
					irow1++;
				}
			} /* for( i = 0; i < voltage_source_count; i++ ) */
		
			if ( irow1 >= 2)
			{
				for( i = 0; i < irow1; i++ )
				{
					isc1 = ipnt[i]-1;
					asmx= geometry.segment_length[isc1];
				
					for (int j = 0; j < neqt; j++ )
						rhs[j] = cplx_00();
				
					rhs[isc1] = cplx_10();
					solves( cm, ip, rhs, neq, 1, geometry.np, geometry.n, geometry.mp, geometry.m, nop, symmetry_array);
					geometry.get_current_coefficients(wavelength, rhs, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
				
					for (int j = 0; j < irow1; j++ )
					{
						isc1= ipnt[j]-1;
						cmn[j+i*ndimn]= rhs[isc1]/ asmx;
					}
				} /* for( i = 0; i < irow1; i++ ) */
			
				asmx=0.0;
				asa=0.0;
			
				for( i = 1; i < irow1; i++ )
				{
					for (int j = 0; j < i; j++ )
					{
						cux = cmn[i+j*ndimn];
						nec_float pwr= abs(( cux- cmn[j+i*ndimn])/ cux);
						asa += pwr* pwr;
				
						if ( pwr >= asmx)
						{
							asmx= pwr;
							nteq= ipnt[i];
							ntsc= ipnt[j];
						}
					} /* for( j = 0; j < i; j++ ) */	
				} /* for( i = 1; i < irow1; i++ ) */
	
				asa= sqrt( asa*2./ (nec_float)( irow1*( irow1-1)));
				fprintf( output_fp, "\n\n"
					"   MAXIMUM RELATIVE ASYMMETRY OF THE DRIVING POINT ADMITTANCE\n"
					"   MATRIX IS %10.3E FOR SEGMENTS %d AND %d\n"
					"   RMS RELATIVE ASYMMETRY IS %10.3E",
					asmx, nteq, ntsc, asa );
			} /* if ( irow1 >= 2) */
		} /* if ( masym != 0) */
	
		/* solution of network equations */
		if ( network_count != 0)
		{
			// zero the cmn array, and the rhnx array
			cmn.fill(cplx_00());
			rhnx.fill(cplx_00());
			
/*			for( i = 0; i < ndimn; i++ )
			{
				rhnx[i]=cplx_00();
				for (int j = 0; j < ndimn; j++ )
					cmn[j+i*ndimn]=cplx_00();
			}
*/			
			nteq=0;
			ntsc=0;
	
			/*	sort network and source data and
				assign equation numbers to segments */
			for (int j = 0; j < network_count; j++ )
			{
				nseg1= iseg1[j];
				nseg2= iseg2[j];
			
				if ( ntyp[j] <= 1)
				{
					y11r= x11r[j];
					y11i= x11i[j];
					y12r= x12r[j];
					y12i= x12i[j];
					y22r= x22r[j];
					y22i= x22i[j];
				}
				else
				{
					y22r= two_pi() * x11i[j]/ wavelength;
					y12r=0.;
					y12i=1./( x11r[j]* sin( y22r));
					y11r= x12r[j];
					y11i=- y12i* cos( y22r);
					y22r= x22r[j];
					y22i= y11i+ x22i[j];
					y11i= y11i+ x12i[j];
				
					if ( ntyp[j] != 2)
					{
						y12r=- y12r;
						y12i=- y12i;
					}
				} /* if ( ntyp[j] <= 1) */
		
				jump1 = false;
				for( i = 0; i < voltage_source_count; i++ )
				{
					if ( nseg1 == source_segment_array[i])
					{
						isc1 = i;
						jump1 = true;
						break;
					}
				}
		
				jump2 = false;
				if ( ! jump1 )
				{
					isc1=-1;
			
					for( i = 0; i < nteq; i++ )
					{
						if ( nseg1 == nteqa[i])
						{
							irow1 = i;
							jump2 = true;
							break;
						}
					}
			
					if ( ! jump2 )
					{
						irow1= nteq;
						nteqa[nteq]= nseg1;
						nteq++;
					}
				} /* if ( ! jump1 ) */
				else
				{
					for( i = 0; i < ntsc; i++ )
					{
						if ( nseg1 == ntsca[i])
						{
							irow1 = ndimn- (i+1);
							jump2 = true;
							break;
						}
					}
			
					if ( ! jump2 )
					{
						irow1= ndimn- (ntsc+1);
						ntsca[ntsc]= nseg1;
						vsrc[ntsc]= source_voltage_array[isc1];
						ntsc++;
					}
			
				} /* if ( ! jump1 ) */
		
				jump1 = false;
				for( i = 0; i < voltage_source_count; i++ )
				{
					if ( nseg2 == source_segment_array[i])
					{
						isc2= i;
						jump1 = true;
						break;
					}
				}
		
				jump2 = false;
				if ( ! jump1 )
				{
					isc2=-1;
			
					for( i = 0; i < nteq; i++ )
					{
						if ( nseg2 == nteqa[i])
						{
							irow2= i;
							jump2 = true;
							break;
						}
					}
			
					if ( ! jump2 )
					{
						irow2= nteq;
						nteqa[nteq]= nseg2;
						nteq++;
					}
				}  /* if ( ! jump1 ) */
				else
				{
					for( i = 0; i < ntsc; i++ )
					{
						if ( nseg2 == ntsca[i])
						{
							irow2 = ndimn- (i+1);
							jump2 = true;
							break;
						}
					}
		
					if ( ! jump2 )
					{
						irow2= ndimn- (ntsc+1);
						ntsca[ntsc]= nseg2;
						vsrc[ntsc]= source_voltage_array[isc2];
						ntsc++;
					}
				} /* if ( ! jump1 ) */
		
				/* fill network equation matrix and right hand side vector with */
				/* network short-circuit admittance matrix coefficients. */
				if ( isc1 == -1)
				{
					cmn[irow1+irow1*ndimn] -= nec_complex( y11r, y11i)* geometry.segment_length[nseg1-1];
					cmn[irow1+irow2*ndimn] -= nec_complex( y12r, y12i)* geometry.segment_length[nseg1-1];
				}
				else
				{
					rhnx[irow1] += nec_complex( y11r, y11i)* source_voltage_array[isc1]/wavelength;
					rhnx[irow2] += nec_complex( y12r, y12i)* source_voltage_array[isc1]/wavelength;
				}
			
				if ( isc2 == -1)
				{
					cmn[irow2+irow2*ndimn] -= nec_complex( y22r, y22i)* geometry.segment_length[nseg2-1];
					cmn[irow2+irow1*ndimn] -= nec_complex( y12r, y12i)* geometry.segment_length[nseg2-1];
				}
				else
				{
					rhnx[irow1] += nec_complex( y12r, y12i)* source_voltage_array[isc2]/wavelength;
					rhnx[irow2] += nec_complex( y22r, y22i)* source_voltage_array[isc2]/wavelength;
				}
			} /* for( j = 0; j < network_count; j++ ) */
	
			/*	add interaction matrix admittance
				elements to network equation matrix */
			for( i = 0; i < nteq; i++ )
			{
				for (int j = 0; j < neqt; j++ )
					rhs[j] = cplx_00();
				
				irow1= nteqa[i]-1;
				rhs[irow1]=cplx_10();
				solves( cm, ip, rhs, neq, 1, geometry.np, geometry.n, geometry.mp, geometry.m, nop, symmetry_array);
				geometry.get_current_coefficients(wavelength, rhs, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
				
				for (int j = 0; j < nteq; j++ )
				{
					irow1= nteqa[j]-1;
					cmn[i+j*ndimn] += rhs[irow1];
				}
			} /* for( i = 0; i < nteq; i++ ) */
		
			/* factor network equation matrix */
			lu_decompose( nteq, cmn, ipnt, ndimn);
			
		} /* if ( network_count != 0) */
	} /* if ( ntsol != 0) */

	if (0 == network_count)
	{
		/* solve for currents when no networks are present */
		solves( cm, ip, einc, neq, 1, geometry.np, geometry.n, geometry.mp, geometry.m, nop, symmetry_array);
		geometry.get_current_coefficients(wavelength, einc, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
		ntsc=0;
	}
	else // if ( network_count != 0)
	{
		/* add to network equation right hand side */
		/* the terms due to element interactions */
		for( i = 0; i < neqt; i++ )
			rhs[i]= einc[i];
	
		solves( cm, ip, rhs, neq, 1, geometry.np, geometry.n, geometry.mp, geometry.m, nop, symmetry_array);
		geometry.get_current_coefficients(wavelength, rhs, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
	
		for( i = 0; i < nteq; i++ )
		{
			irow1= nteqa[i]-1;
			rhnt[i]= rhnx[i]+ rhs[irow1];
		}

		/* solve network equations */
		solve( nteq, cmn, ipnt, rhnt, ndimn);
	
		/* add fields due to network voltages to electric fields */
		/* applied to structure and solve for induced current */
		for( i = 0; i < nteq; i++ )
		{
			irow1= nteqa[i]-1;
			einc[irow1] -= rhnt[i];
		}
	
		solves( cm, ip, einc, neq, 1, geometry.np, geometry.n, geometry.mp, geometry.m, nop, symmetry_array);
		geometry.get_current_coefficients(wavelength, einc, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);

		if ( nprint == 0)
		{
			fprintf( output_fp, "\n\n\n"
				"                          "
				"--------- STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS --------" );
		
			fprintf( output_fp, "\n"
				"  TAG   SEG       VOLTAGE (VOLTS)          CURRENT (AMPS)        "
				" IMPEDANCE (OHMS)       ADMITTANCE (MHOS)     POWER\n"
				"  No:   No:     REAL      IMAGINARY     REAL      IMAGINARY    "
				" REAL      IMAGINARY     REAL      IMAGINARY   (WATTS)" );
		}

		for( i = 0; i < nteq; i++ )
		{
			int segment_number = nteqa[i];
			int segment_index = segment_number-1;
			nec_complex voltage = rhnt[i]* geometry.segment_length[segment_index]* wavelength;
			nec_complex current = einc[segment_index]* wavelength;
			nec_complex admittance = current / voltage;
			nec_complex impedance = voltage / current;
			int segment_tag = geometry.segment_tags[irow1];
			nec_float power = em::power(voltage,current);
			network_power_loss= network_power_loss - power;
			
			if ( nprint == 0)
				fprintf( output_fp, "\n"
					" %4d %5d %11.4E %11.4E %11.4E %11.4E"
					" %11.4E %11.4E %11.4E %11.4E %11.4E",
					segment_tag, segment_number, real(voltage), imag(voltage), real(current), imag(current),
					real(impedance), imag(impedance), real(admittance), imag(admittance), power );
		}

		for( i = 0; i < ntsc; i++ )
		{
			irow1= ntsca[i]-1;
			vlt= vsrc[i];
			cux= einc[irow1]* wavelength;
			ymit= cux/ vlt;
			zped= vlt/ cux;
			irow2= geometry.segment_tags[irow1];
			
			nec_float pwr= em::power(vlt,cux);
			network_power_loss= network_power_loss- pwr;
		
			if ( nprint == 0)
				fprintf( output_fp, "\n"
					" %4d %5d %11.4E %11.4E %11.4E %11.4E"
					" %11.4E %11.4E %11.4E %11.4E %11.4E",
					irow2, irow1+1, real(vlt), imag(vlt), real(cux), imag(cux),
					real(zped), imag(zped), real(ymit), imag(ymit), pwr );
		} /* for( i = 0; i < ntsc; i++ ) */
	} /* if ( network_count != 0) */

	if ( (voltage_source_count+nvqd) == 0)
		return;
	
	nec_antenna_input* antenna_input = new nec_antenna_input();
	s_results.add(antenna_input);
	
	s_output.end_section();
	fprintf( output_fp, 
		"                        "
		"--------- ANTENNA INPUT PARAMETERS ---------" );
	
	fprintf( output_fp, "\n"
		"  TAG   SEG       VOLTAGE (VOLTS)         "
		"CURRENT (AMPS)         IMPEDANCE (OHMS)    "
		"    ADMITTANCE (MHOS)     POWER\n"
		"  NO.   NO.     REAL      IMAGINARY"
		"     REAL      IMAGINARY     REAL      "
		"IMAGINARY    REAL       IMAGINARY   (WATTS)" );
	
	for( i = 0; i < voltage_source_count; i++ )
	{
		int segment_index = source_segment_array[i]-1;
		nec_complex voltage = source_voltage_array[i];
		nec_complex current = einc[segment_index] * wavelength;
		
		bool add_as_network_loss = false;
		
		// the following loop is completely mysterious!
		for (int j = 0; j < ntsc; j++ )
		{
			// I am now almost sure that the following code is not correct.
			// This modifies the current, however if the inner loop is executed more
			// than once, then only the last current modification is kept!
			 
			if ( ntsca[j] == segment_index+1)
			{
				int row_index = ndimn - (j+1);
				int row_offset = row_index*ndimn;
				
				// I wish I knew what was going on here...
				nec_complex temp = rhnx[row_index]; // renamed current -> temp to avoid confusion
				for (int k = 0; k < nteq; k++ )
					temp -= cmn[k + row_offset]*rhnt[k];
					
				current = (temp + einc[segment_index])* wavelength;
				add_as_network_loss = true;
					
#warning "This loop is messed up. The j is inside another j loop"
				// I have removed the j from the "for (int k = 0; k < nteq; k++ )" loop 
				// and placed this"j=nteq" statement here.
				j = nteq;
			}
		}
			
		nec_complex admittance = current / voltage;
		nec_complex impedance = voltage / current;
		nec_float power = em::power(voltage,current);
		
		if ( add_as_network_loss )
			network_power_loss += power;
			
		input_power += power;
		
		int segment_tag = geometry.segment_tags[segment_index];
	
		antenna_input->set_input(
			segment_tag, segment_index+1,
			voltage, current, impedance, admittance, power);
		
		fprintf( output_fp,	"\n"
			" %4d %5d %11.4E %11.4E %11.4E %11.4E"
			" %11.4E %11.4E %11.4E %11.4E %11.4E",
			segment_tag, segment_index+1, real(voltage), imag(voltage), real(current), imag(current),
			real(impedance), imag(impedance), real(admittance), imag(admittance), power );
		
	} /* for( i = 0; i < voltage_source_count; i++ ) */

	
	for( i = 0; i < nvqd; i++ )
	{
		int segment_index = ivqd[i]-1;
		nec_complex voltage = vqd[i];
		
		nec_complex _ai( air[segment_index], aii[segment_index]);
		nec_complex _bi( bir[segment_index], bii[segment_index]);
		nec_complex _ci( cir[segment_index], cii[segment_index]);
		
		// segment length is measured in wavelengths. The pase is therefore the length in wavelengths
		// multiplied by pi().
		nec_float segment_length_phase = geometry.segment_length[segment_index] * pi(); // TCAM CHANGED TO pi() (from TP*.5)!!
		
		nec_complex current = ( _ai - _bi* sin(segment_length_phase)+ _ci * cos(segment_length_phase)) * wavelength;
		
		nec_complex admittance = current / voltage;
		nec_complex impedance = voltage / current;
		nec_float power = em::power(voltage,current);
		
		input_power += power;
		
		int segment_tag = geometry.segment_tags[segment_index];
	
		antenna_input->set_input(
			segment_tag, segment_index+1,
			voltage, current, impedance, admittance, power);
		
		fprintf( output_fp,	"\n"
			" %4d %5d %11.4E %11.4E %11.4E %11.4E"
			" %11.4E %11.4E %11.4E %11.4E %11.4E",
			segment_tag, segment_index+1, real(voltage), imag(voltage), real(current), imag(current),
			real(impedance), imag(impedance), real(admittance), imag(admittance), power );
		
	} /* for( i = 0; i < nvqd; i++ ) */
}
