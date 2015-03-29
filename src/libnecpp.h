/***************************************************************************
 *   Copyright (C) 2004-2008,2015 by Tim Molteno                                     *
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
/**
 * \file libnecpp.h
 * \brief nec++ Library Functions.
 * \section _how_to_use How to use libNEC. 
 * Have a look at \ref test_nec.c
 * 
 * \example test_nec.c
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

/*! \brief Create an nec_context and initialize it.

\par Note: Do NOT delete or free the nec_context yourself, rather call nec_delete() to free memory associated with the nec simulation.
*/
nec_context* nec_create(void);

/*!\brief Delete an nec_context object. 
 */
long nec_delete(nec_context* in_context);

/*!\brief Benchmark the libnecpp engine. A score of 1 is roughly an Athlon XP 1800. 
 */
long nec_benchmark(void);


/*! \brief Generates segment geometry for a straigt wire
    \param in_context The nec_context created with nec_create()
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
    
    \remark All co-ordinates are in meters.
*/
long nec_wire(nec_context* in_context, int tag_id, int segment_count,
              double xw1, double yw1, double zw1,
              double xw2, double yw2, double zw2, 
              double rad, double rdel, double rrad);


/*! \brief Surface Patch (SP Card)
    \param in_context The nec_context created with nec_create()
    \param ns The Patch Type.
          \arg \c 0 (default) arbitrary patch shape
          \arg \c 1 rectangular patch
          \arg \c 2 triangular patch
          \arg \c 3 quadrilateral patch 
    \param x1 The x coordinate of patch corner1.
    \param y1 The y coordinate of patch corner1.
    \param z1 The z coordinate of patch corner1.
    \param x2 The x coordinate of patch corner2.
    \param y2 The y coordinate of patch corner2.
    \param z2 The z coordinate of patch corner2.
    
    \remark All co-ordinates are in meters, except for arbitrary patches where the angles are in degrees
*/
long nec_sp_card(nec_context* in_context, int ns,
    double x1, double y1, double z1,
    double x2, double y2, double z2);

/*! \brief Surface Patch Continuation (SC Card)
    \param in_context The nec_context created with nec_create()
    \param i2  Weird integer parameter.
    \param x3 The x coordinate of patch corner 3.
    \param y3 The y coordinate of patch corner 3.
    \param z3 The z coordinate of patch corner 3.
    \param x4 The x coordinate of patch corner 4.
    \param y4 The y coordinate of patch corner 4.
    \param z4 The z coordinate of patch corner 4.
    
    \remark All co-ordinates are in meters.
*/
long nec_sc_card(nec_context* in_context, int i2,
    double x3, double y3, double z3,
    double x4, double y4, double z4);

/*! \brief Coordinate Transformation
 * 
 * \param itsi  Tag number increment.
 * \param nprt  The number of new Structures to be generated
 * \param ROX   Angle in degrees through which the structure is rotated about
 *              the X-axis.  A positive angle causes a right-hand rotation.
 * \param ROY   Angle of rotation about Y-axis.
 * \param ROZ   Angle of rotation about
 * \param XS    X, Y. Z components of vector by which
 * \param YS    structure is translated with respect to
 * \param ZS    the coordinate system.
 * \param ITS   This number is input as a decimal number but is rounded
 *             to an integer before use.  Tag numbers are searched sequentially
 *             until a segment having a tag of this segment through the end of
 *             the sequence of segments is moved by the card.  If ITS is zero 
 *             the entire structure is moved.
**/
long nec_gm_card(nec_context* in_context, int itsi, int nrpt,
                 double rox, double roy, double roz, double xs,
                 double ys, double zs, int its );

/*!\brief Parameters:
        Integers
   \param i1 - Tag number increment.
   \param i2 - This integer is divided into three independent digits, in
                  columns 8, 9, and 10 of the card, which control reflection
                  in the three orthogonal coordinate planes.  A one in column
                  8 causes reflection along the X-axis (reflection in Y, Z
                  plane); a one in column 9 causes reflection along the Y-axis;
                  and a one in column 10 causes reflection along the Z axis.
                  A zero or blank in any of these columns causes the corres-
                  ponding reflection to be skipped.

   \remark Any combination of reflections along the X, Y and Z axes may be used. 
   For example, 101 for (I2) will cause reflection along axes X and Z, and 111 will 
   cause reflection along axes X, Y and Z. When combinations of reflections are requested, 
   the reflections are done in reverse alphabetical order. That is, if a structure is 
   generated in a single octant of space and a GX card is then read with I2 equal to 111, 
   the structure is first reflected along the Z-axis; the structure and its image are 
   then reflected along the Y-axis; and, finally, these four structures are reflected 
   along the X-axis to fill all octants. This order determines the position of a segment 
   in the sequence and, hence, the absolute segment numbers.
   
   \remark The tag increment I1 is used to avoid duplication of tag numbers in the image 
   segments. All valid tags on the original structure are incremented by I1 on the image.
   When combinations of reflections are employed, the tag increment is doubled after each 
   reflection. Thus, a tag increment greater than or equal to the largest tag an the 
   original structure will ensure that no duplicate tags are generated. For example, 
   if tags from 1 to 100 are used on the original structure with I2 equal to 011 and 
   a tag increment of 100, the first reflection, along the Z-axis, will produce tags 
   from 101 to 200; and the second reflection, along the Y-axis, will produce tags 
   from 201 to 400, as a result of the increment being doubled to 200. 
 */
long nec_gx_card(nec_context* in_context, int i1, int i2);


/*! \brief Indicate that the geometry is complete (GE card)
 * \param in_context The nec_context created with nec_create()
 * \param gpflag Geometry ground plain flag.
 *    \arg \c 0 - no ground plane is present.
 *    \arg \c 1 - Indicates a ground plane is present. Structure symmetry is modified as required, and the current expansion is modified so that the currents an segments touching the ground (x, Y plane) are interpolated to their images below the ground (charge at base is zero)
 *    \arg \c -1 - indicates a ground is present. Structure symmetry is modified as required. Current expansion, however, is not modified, Thus, currents on segments touching the ground will go to zero at the ground. 
 * \param card_int_2 Unused (set to zero)
 **/
long nec_geometry_complete(nec_context* in_context, int gpflag, int card_int_2);


/*! \brief Get the last error message
 * All functions return a long. If this is != 0. Then an error has occurred.
 * The error message can be retrieved with this function.
 **/
const char* nec_error_message(void);

/*
  NEC card functions.
*/
/*!\brief Ground Card
  Examples:

  1) Infinite ground plane
    nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);

  2) Radial Wire Ground Plane (4 wires, 2 meters long, 5mm in radius)
    nec_gn_card(nec, 4, 0, 0.0, 0.0, 2.0, 0.005, 0.0, 0.0)
    
  \param iperf Ground-type flag
  \arg \c -1 Nullifies ground parameters previously used and sets free-space condition. The remainder of the parameters are ignored in this case.
  \arg \c 0 Finite ground, reflection coefficient approximation
  \arg \c 1 Perfectly conducting ground.
  \arg \c 2 Finite ground, Sommerfeld/Norton method.
  
  \param nradl Number of radial wires in the ground screen approximation, O implies no ground screen.
  
  \param epse Relative dielectric constant for ground in the vicinity of the antenna. Zero in the case of perfect ground.
  \param sig Conductivity in mhos/meter of the ground in the vicinity of the antenna. Use zero in the case of a perfect ground. If SIG is input as a negative number, the complex dielectric constant Ec = Er -j sigma/omaga epslon is set to EPSR - |SIG|. 
*/
long nec_gn_card(nec_context* in_context, int iperf, int nradl, double epse, double sig, double tmp3, double tmp4, double tmp5, double tmp6);

/*! \brief FR card
 * \param in_context The nec_context created with nec_create()
 * \param in_ifrq 0 is a linear range of frequencies, 1 is a log range.
 * \param in_nfrq The number of frequencies
 * \param in_freq_mhz The starting frequency in MHz.
 * \param in_del_freq The frequency step (in MHz for ifrq = 0)
 */
long nec_fr_card(nec_context* in_context, int in_ifrq, int in_nfrq, double in_freq_mhz, double in_del_freq);



/*! \brief LD card (Loading)
* \param in_context The nec_context created with nec_create()
* \param ldtyp Type of loading (5 = segment conductivity)
* \param ldtag Tag (zero for absolute segment numbers, or in conjunction with 0 for next parameter, for all segments)
* \param ldtagf Equal to m specifies the mth segment of the set of segments whose tag numbers equal the tag number 
* specified in the previous parameter. If the previous parameter (LDTAG) is zero, LDTAGF then specifies an absolute segment number. 
* If both LDTAG and LDTAGF are zero, all segments will be loaded. 
* \param ldtagt Equal to n specifies the nth segment of the set of segments whose tag numbers equal the tag number specified 
* in the parameter LDTAG. This parameter must be greater than or equal to the previous parameter. 
* The loading specified is applied to each of the mth through nth segments of the set of segments having tags 
* equal to LDTAG. Again if LDTAG is zero, these parameters refer to absolute segment numbers. 
* If LDTAGT is left blank, it is set equal to the previous parameter (LDTAGF).
* \remark Floating Point Input for the Various Load Types:
*/
long nec_ld_card(nec_context* in_context, int ldtyp, int ldtag, int ldtagf, int ldtagt, double tmp1, double tmp2, double tmp3);


long nec_ex_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
long nec_tl_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
long nec_nt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);

/*! \brief XQ Card (Execute)
 * 
 * Purpose:   To cause program execution at points in the data stream where
 *            execution is not automatic.  Options on the card also allow for
 *            automatic generation of radiation patterns in either of two vertical
 *            cuts.
 * \param in_context The nec_context created with nec_create()
 * \param itmp1 Options controlled by (I1) are:
 *            0 - no patterns requested (normal case).
 *            1 - generates a pattern cut in the XZ plane, i.e., phi = 0 degrees
 *                and theta varies from 0 degrees to 90 degrees in 1 degree steps.
 *            2 - generates a pattern cut in the YZ plane, i.e., phi = 90 degrees
 *                theta varies from 0 degrees to 90 degrees in 1 degree steps.
 *            3 - generates both of the cuts described for the values 1 and 2.
 **/
long nec_xq_card(nec_context* in_context, int itmp1);

long nec_gd_card(nec_context* in_context, double tmp1, double tmp2, double tmp3, double tmp4);

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
long nec_rp_card(nec_context* in_context,
	int calc_mode, int n_theta, int n_phi,
	int output_format, int normalization, int D, int A,	
	double theta0, double phi0, double delta_theta, double delta_phi,
	double radial_distance, double gain_norm);

/*! \brief Print Flag (Printing of Currents
 * \param IPTFLG Print control flag, specifies the type of format used in printing segment currents. The options are:
      \arg \c -2 - all currents printed. This it a default value for the program if the card is Omitted.
      \arg \c -1 - suppress printing of all wire segment currents.
      \arg \c O - current printing will be limited to the segments specified by the next three parameters.
      \arg \c 1 - currents are printed by using a format designed for a receiving pattern (refer to output section in this manual Only currents for the segments specified by the next three parameters are printed.
      \arg \c 2 - same as for 1 above; in addition, however, the current for one Segment will Cue normalized to its maximum, ant the normalized values along with the relative strength in tB will be printed in a table. If the currents for more than one segment are being printed, only currents from the last segment in the group appear in the normalized table.
      \arg \c 3 - only normalized currents from one segment are printed for the receiving pattern case. 

    \param IPTAG - Tag number of the segments for which currents will be printed. 

    \param IPTAGF - Equal to m, specifies the mth segment of the set of segments having the tag numbers of IPTAG, at which printing of currents starts. If IPTAG is zero or blank, then IPTAGF refers to an absolute segment number. If IPTAGF is blank, the current is printed for all segments.

    \param IPTAGT - Equal to n specifies the nth segment of the set of segments having tag numbers of IPTAG. Currents are printed for segments having tag number IPTAG starting at the m th segment in the set and ending at the nth segment. If IPTAG is zero or blank, then IPTAGF and IPTAGT refer to absoulte segment numbers. In IPTAGT is left blank, it is set to IPTAGF.
 */
long nec_pt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);

long nec_pq_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);

long nec_kh_card(nec_context* in_context, double tmp1);
long nec_ne_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
long nec_nh_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);

/*!\brief To control use of the extended thin-wire kernel approximation.
 * \param itmp1 
 * \arg \c -1 Return to normal kernel
 * \arg \c 0 Use Extended thin wire kernel
 * */
long nec_ek_card(nec_context* in_context, int itmp1);

long nec_cp_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
long nec_pl_card(nec_context* in_context, char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4);


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

/*! \brief Impedance: Real Part 
 */
double nec_impedance_real(nec_context* in_context, int freq_index);

/*! \brief Impedance: Imaginary Part 
 */
double nec_impedance_imag(nec_context* in_context, int freq_index);



#ifdef __cplusplus
}
#endif


#endif /* __libnecpp__ */
