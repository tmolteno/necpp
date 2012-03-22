/***************************************************************************
 *   Copyright (C) 2004-2008 by Tim Molteno                                     *
 *   tim@molteno.net                                               *
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
/*! \file libnecpp.h
    \brief nec++ Library Functions.
    
\verbatim
How to use libNEC. 
	
	Enter the following file into test_nec.c, and compile with
	gcc -o test_nec test_nec.c -lnecpp
	
	#include "libnecpp.h"
	#include <stdio.h>
	
	int main(int argc, char **argv)
	{
		nec_context* nec;
		double gain;
			
		nec = nec_create();
		nec_wire(nec, 0, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 21, -0.042, 0.008, 0.017, -0.048, 0.021, -0.005, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 70, -0.048, 0.021, -0.005, 0.039, 0.032, -0.017, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 70, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 50, -0.042, 0.008, 0.017, 0.017, -0.015, 0.014, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 85, -0.027, 0.04, -0.031, 0.046, -0.01, 0.028, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 70, 0.017, -0.015, 0.014, -0.048, -0.038, -0.04, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 77, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0);
		nec_geometry_complete(nec, 0, 0);
		
		nec_gn_card(nec, -1,0,0.0, 0.0, 0.0,0.0, 0.0, 0.0);
		nec_ld_card(nec, 5,0,0,0,3.72e7,0.0,0.0);
		nec_pt_card(nec, -1, 0, 0, 0);
		nec_ex_card(nec, 1, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		nec_fr_card(nec, 0, 2, 2400.0, 100.0);
		nec_rp_card(nec, 0, 1, 1, 0,5,0,0, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0);

		printf("Impedance: %f, %f\n",nec_impedance_real(nec,0), nec_impedance_imag(nec,0));
		printf("Gain: %f, %f +/- %f dB\n",nec_gain_max(nec,0), nec_gain_mean(nec,0), nec_gain_sd(nec,0));
		printf("RHCP Gain: %f, %f +/- %f dB\n",nec_gain_rhcp_max(nec,0), nec_gain_rhcp_mean(nec,0), nec_gain_rhcp_sd(nec,0));
		printf("LHCP Gain: %f, %f +/- %f dB\n",nec_gain_lhcp_max(nec,0), nec_gain_lhcp_mean(nec,0), nec_gain_lhcp_sd(nec,0));
		
		nec_delete(nec);
		
		
		return 0;
	}
\endverbatim
*/

#ifndef __libnecpp__
#define __libnecpp__

/*! A Struct to represent the nec_context class.*/
#ifdef __cplusplus
class nec_context;
#else
struct nec_context;
typedef struct nec_context nec_context;
#endif


#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Construct and initialize an nec_context */
nec_context* nec_create();

/*!\brief Delete an nec_context object. */
long nec_delete(nec_context* in_context);

/*!\brief Benchmark the libnecpp engine. A score of 100 is roughly an Athlon XP 1800. */
long nec_benchmark();


/*! \brief Generates segment geometry for a straigt wire
	\param in_context The nec_context created with nec_create()
	\param tag_id
	\param segment_count Number of Elements (should be around 12-20 per wavelength)
	\param rad Wire radius of first segment (in Meters)
	\param rdel Ratio of the length of a segment to the length of the previous segment.  (Set to 1.0 if segments have uniform length)
	\param rrad The ratio of the radii of adjacent segments (Set to 1.0 if not tapered)
*/
/*! Add a wire to the geometry,

All co-ordinates are in meters.

	\param tag_id The tag ID.
	\param segment_count The number of segments.
	\param xw1 The x coordinate of the wire starting point.
	\param yw1 The y coordinate of the wire starting point.
	\param zw1 The z coordinate of the wire starting point.
	\param xw2 The x coordinate of the wire ending point.
	\param yw2 The y coordinate of the wire ending point.
	\param zw2 The z coordinate of the wire ending point.
	\param rad The wire radius (meters)
	\param rdel For tapered wires, the. Otherwise set to 1.0
	\param rrad For tapered wires, the. Otherwise set to 1.0
*/
void nec_wire(nec_context* in_context, int tag_id, int segment_count,
		double xw1, double yw1, double zw1,
		double xw2, double yw2, double zw2, 
		double rad, double rdel, double rrad);

/*! \brief Indicate that the geometry is complete (GE card)
	\param in_context The nec_context created with nec_create()
	\param gpflag Geometry ground plain flag.
		0 - no ground plane is present. 
		1 - Indicates a ground plane is present. Structure symmetry is modified as required, and the current expansion is modified so that the currents an segments touching the ground (x, Y plane) are interpolated to their images below the ground (charge at base is zero) 
		-1 - indicates a ground is present. Structure symmetry is modified as required. Current expansion, however, is not modified, Thus, currents on segments touching the ground will go to zero at the ground. 
	\param card_int_2 Unused (set to zero)
*/
void nec_geometry_complete(nec_context* in_context, int gpflag, int card_int_2);


/*
	NEC card functions.
*/
/*!\brief Ground Card
Examples:

1) Infinite ground plane
	nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);

2) Radial Wire Ground Plane (4 wires, 2 meters long, 5mm in radius)
	nec_gn_card(nec, 4, 0, 0.0, 0.0, 2.0, 0.005, 0.0, 0.0)
*/
void nec_gn_card(nec_context* in_context, int itmp1, int itmp2, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);

/*!
 * FR card
 *	@param in_context The nec_context created with nec_create()
 *	@param in_ifrq 0 is a linear range of frequencies, 1 is a log range.
 *	@param in_nfrq The number of frequencies
 *	@param in_freq_mhz The starting frequency in MHz.
 *	@param in_del_freq The frequency step (in MHz for ifrq = 0)
 */
void nec_fr_card(nec_context* in_context, int in_ifrq, int in_nfrq, double in_freq_mhz, double in_del_freq);



/*!
* LD card (Loading)
*	@param in_context The nec_context created with nec_create()
*	@param ldtyp Type of loading (5 = segment conductivity)
*	@param ldtag Tag (zero for absolute segment numbers, or in conjunction with 0 for next parameter, for all segments)
*	@param ldtagf Equal to m specifies the mth segment of the set of segments whose tag numbers equal the tag number specified in the previous parameter. If the previous parameter (LDTAG) is zero, LDTAGF then specifies an absolute segment number. If both LDTAG and LDTAGF are zero, all segments will be loaded. 
*	@param ldtagt Equal to n specifies the nth segment of the set of segments whose tag numbers equal the tag number specified in the parameter LDTAG. This parameter must be greater than or equal to the previous param- eter. The loading specified is applied to each of the mth through nth segments of the set of segments having tags equal to LDTAG. Again if LDTAG is zero, these parameters refer to absolute segment numbers. If LDTAGT is left blank, it is set equal to the previous parameter (LDTAGF).

Floating Point Input for the Various Load Types:
*/
void nec_ld_card(nec_context* in_context, int ldtyp, int ldtag, int ldtagf, int ldtagt, double tmp1, double tmp2, double tmp3);


void nec_ex_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_tl_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_nt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_xq_card(nec_context* in_context, int itmp1);
void nec_gd_card(nec_context* in_context, double tmp1, double tmp2, double tmp3, double tmp4);

	/*! \brief Standard radiation pattern parameters 
	
	\param calc_mode This integer selects the mode of calculation for the radiated field. Some values of (calc_mode) will affect the meaning of the remaining parameters on the card. Options available for calc_mode are:
		\arg \c O - normal mode. Space-wave fields are computed. An infinite ground plane is included if it has been specified previously on a GN card; otherwise, the antenna is in free space.
		\arg \c 1 - surface wave propagating along ground is added to the normal space wave. This option changes the meaning of some of the other parameters on the RP card as explained below, and the results appear in a special output format. Ground parameters must have been input on a GN card. The following options cause calculation of only the space wave but with special ground conditions. Ground conditions include a two-medium ground (cliff where the media join in a circle or a line), and a radial wire ground screen. Ground parameters and dimensions must be input on a GN or GD card before the RP card is read. The RP card only selects the option for inclusion in the field calculation. (Refer to the GN and GD cards for further explanation.)
		\arg \c  2 - linear cliff with antenna above upper level. Lower medium parameters are as specified for the second medium on the GN card or on the GD card.
		\arg \c 3 - circular cliff centered at origin of coordinate system: with antenna above upper level. Lower medium parameters are as specified for the second medium on the GN card or on the GD card.
		\arg \c 4 - radial wire ground screen centered at origin.
		\arg \c 5 - both radial wire ground screen and linear cliff.
		\arg \c 6 - both radial wire ground screen ant circular cliff.
	
	\param n_theta The number of theta angles.
	\param n_phi The number of phi angles. 
		 

		 
\param output_format The output format:
	\arg \c 0 major axis, minor axis and total gain printed. 
	\arg \c 1 vertical, horizontal ant total gain printed.
	
\param normalization Controls the type of normalization of the radiation pattern
	\arg \c 0 no normalized gain. 
	\arg \c 1 major axis gain normalized. 
	\arg \c 2 minor axis gain normalized. 
	\arg \c 3 vertical axis gain normalized. 
	\arg \c 4 horizontal axis gain normalized. 
	\arg \c 5 total gain normalized.

\param D Selects either power gain or directive gain for both standard printing and normalization. If the structure excitation is an incident plane wave, the quantities printed under the heading "gain" will actually be the scattering cross section (a/lambda 2 ) and will not be affected by the value of d. The column heading for the output will still read "power" or "directive gain," however. 
	\arg \c 0 power gain. 
	\arg \c 1 directive gain.


\param A - Requests calculation of average power gain over the region covered by field points. 
	\arg \c 0 no averaging. 
	\arg \c 1 average gain computed. 
	\arg \c 2 average gain computed, printing of gain at the field points used for averaging is suppressed. If n_theta or NPH is equal to one, average gain will not be computed for any value of A since the area of the region covered by field points vanishes.


	
	\param theta0 - Initial theta angle in degrees (initial z coordinate in meters if calc_mode = 1).
	
	\param phi0 - Initial phi angle in degrees.
	
	\param delta_theta - Increment for theta in degrees (increment for z in meters if calc_mode = 1).
	
	\param delta_phi - Increment for phi in degrees.
	
	\param radial_distance - Radial distance (R) of field point from the origin in meters. radial_distance is optional. If it is zero, the radiated electric field will have the factor exp(-jkR)/R omitted. If a value of R is specified, it should represent a point in the far-field region since near components of the field cannot be obtained with an RP card. (If calc_mode = 1, then radial_distance represents the cylindrical coordinate phi in meters and is not optional. It must be greater than about one wavelength.)
	
	\param gain_norm - Determines the gain normalization factor if normalization has been requested in the normalization parameter. If gain_norm is zero, the gain will be normalized to its maximum value. If gain_norm is not zero, the gain wi11 be normalized to the value of gain_norm.
	
	\remark
	The field point is specified in spherical coordinates (R, sigma, theta), except when the surface wave is computed. For computing the surface wave field (calc_mode = l), cylindrical coordinates (phi, theta, z) are used to accurately define points near the ground plane at large radial distances.
		 
	\remark
	The rp_card() function allows automatic stepping of the field point to compute the field over a region about the antenna at uniformly spaced points.
	\remark
	The integers n_theta and n_phi and floating point numbers theta0, phi0, delta_theta, delta_phi, radial_distance, and gain_norm control the field-point stepping.
		
	\li The nec_rp_card() function will cause the interaction matrix to be computed and factored and the structure currents to be computed if these operations have not already been performed. Hence, all required input parameters must be set before the nec_rp_card() function is called. 
	\li At a single frequency, any number of nec_rp_card() calls may occur in sequence so that different field-point spacings may be used over different regions of space. If automatic frequency stepping is being used (i.e., in_nfrq on the nec_fr_card() function is greater than one), only one nec_rp_card() function will act as data inside the loop. Subsequent calls to nec_rp_card() will calculate patterns at the final frequency. 
	\li When both n_theta and n_phi are greater than one, the angle theta (or Z) will be stepped faster than phi. 
	\li When a ground plane has been specified, field points should not be requested below the ground (theta greater than 90 degrees or Z less than zero.)
	
	*/
void nec_rp_card(nec_context* in_context,
	int calc_mode, int n_theta, int n_phi,
	int output_format, int normalization, int D, int A,	
	double theta0, double phi0, double delta_theta, double delta_phi,
	double radial_distance, double gain_norm);

void nec_pt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
void nec_pq_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
void nec_kh_card(nec_context* in_context, double tmp1);
void nec_ne_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_nh_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_ek_card(nec_context* in_context, int itmp1);
void nec_cp_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
void nec_pl_card(nec_context* in_context, char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4);


/*! \brief Get statistics of the gains in dB.

This function requires a previous rp_card() method to have been called (with the gain normalization set to 5)
\param index The rp_card frequency index
\return The maximum gain in dB or -999.0 if no radiation pattern had been previously requested.
*/
double nec_gain_max(nec_context* in_context, int freq_index);
double nec_gain_min(nec_context* in_context, int freq_index);
double nec_gain_mean(nec_context* in_context, int freq_index);
double nec_gain_sd(nec_context* in_context, int freq_index);

double nec_gain_rhcp_max(nec_context* in_context, int freq_index);
double nec_gain_rhcp_min(nec_context* in_context, int freq_index);
double nec_gain_rhcp_mean(nec_context* in_context, int freq_index);
double nec_gain_rhcp_sd(nec_context* in_context, int freq_index);

double nec_gain_lhcp_max(nec_context* in_context, int freq_index);
double nec_gain_lhcp_min(nec_context* in_context, int freq_index);
double nec_gain_lhcp_mean(nec_context* in_context, int freq_index);
double nec_gain_lhcp_sd(nec_context* in_context, int freq_index);

double nec_impedance_real(nec_context* in_context, int freq_index);
double nec_impedance_imag(nec_context* in_context, int freq_index);



#ifdef __cplusplus
}
#endif


#endif /* __libnecpp__ */
