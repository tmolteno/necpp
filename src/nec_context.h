#ifndef __nec_context__
#define __nec_context__

/*
  Copyright (C) 2004-2008,2015  Timothy C.A. Molteno
  tim@molteno.net 
  
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

#include "common.h"
#include "c_ggrid.h"
#include "math_util.h"
#include "matrix_algebra.h"
#include "electromag.h"
#include "nec_radiation_pattern.h"
#include "nec_results.h"
#include "nec_structure_currents.h"
#include "nec_output.h"
#include "nec_ground.h"
#include "c_plot_card.h"

class c_geometry;

enum excitation_return
{
    FREQ_PRINT_NORMALIZATION = 0,
    FREQ_LOOP_CONTINUE = 1,
    FREQ_LOOP_CARD_CONTINUE = 2
};

/*! \brief Using nec_context
 * \file nec_context.h
 * 
 * The following code shows an example of how the nec_context class is used.
 * 
 * \example test_cpp.cpp
 */


/*! \brief Container for an nec2++ simulation
 * 
 * An nec_context object is the container for an nec2++ simulation. A c_geometry object
 * is associated with the nec_context, and then after the simulation is done, the results
 * can be requested from this object.
 * */

class nec_context
{
public:
  nec_context();
  virtual ~nec_context();
  
  // Called after construction...
  void initialize();
  

  void calc_prepare();
  
  inline c_geometry* get_geometry()  {
      return m_geometry;
  }
  
  /*! \brief Get the maximum gain in dB.
  
  This function requires a previous rp_card() method to have been called (with gain normalization requested)
  
  \return The maximum gain in dB or -999.0 if no radiation pattern had been previously requested.
  */
  double get_gain(int freq_index, int theta_index, int phi_index)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)       return -999.0;  
          return rp->get_power_gain(theta_index, phi_index);
  }

  double get_gain_max(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)       return -999.0;  
          return rp->get_gain_max();
  }

  double get_gain_min(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_min();
  }

  double get_gain_mean(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_mean();
  }
  
  double get_gain_sd(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_sd();
  }
  
  /********************** RHCP ********************************/
  double get_gain_rhcp_max(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_rhcp_max();
  }

  double get_gain_rhcp_min(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_rhcp_min();
  }

  double get_gain_rhcp_mean(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_rhcp_mean();
  }
  
  double get_gain_rhcp_sd(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_rhcp_sd();
  }
  
  /********************** LHCP ********************************/
  double get_gain_lhcp_max(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_lhcp_max();
  }

  double get_gain_lhcp_min(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_lhcp_min();
  }

  double get_gain_lhcp_mean(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_lhcp_mean();
  }
  
  double get_gain_lhcp_sd(int freq_index = 0)  {
          nec_radiation_pattern* rp = get_radiation_pattern(freq_index);
          if (NULL == rp)  return -999.0;  
          return rp->get_gain_lhcp_sd();
  }
  
  /****************** IMPEDANCE CHARACTERISTICS *********************/
  
  /*! \brief Impedance: Real Part */
  double get_impedance_real(int freq_index = 0)  {
          nec_antenna_input* ipt = get_input_parameters(freq_index);
          if (NULL == ipt) return -999.0;
          vector<nec_complex>& imp(ipt->get_impedance());
          return imp.back().real();
  }
  /*! \brief Impedance: Imaginary Part */
  double get_impedance_imag(int freq_index = 0)  {
          nec_antenna_input* ipt = get_input_parameters(freq_index);
          if (NULL == ipt) return -999.0;  
          vector<nec_complex>& imp(ipt->get_impedance());
          return imp.back().imag();
  }
          
  /*! \brief Get Antenna Input Parameter Results
  \param index The zero-based index for the result (simulations can return more than one set of results).
  \return The requested antenna input parameter data (or NULL if the result does not exist).
  \note You must NOT delete the nec_antenna_input object when finished with it.
  */
  inline nec_antenna_input* get_input_parameters(int index)  {
          return m_results.get_antenna_input(index);
  }
  
  /*! \brief Get Normalized Receiving Pattern Results
  \param index The zero-based index for the result (simulations can return more than one set of results).
  \return The requested radiation pattern data (or NULL if the result does not exist).
  \note You must NOT delete the nec_norm_rx_pattern object when finished with it.
  */
  inline nec_norm_rx_pattern* get_norm_rx_pattern(int index)  {
          return m_results.get_norm_rx_pattern(index);
  }
  
  /*! \brief Get Radiation Pattern results
  \param index The zero-based index for the result (simulations can return more than one set of results).
  \return The requested radiation pattern data (or NULL if the result does not exist).
  \note You must NOT delete the results object when finished with it.
  */
  inline nec_radiation_pattern* get_radiation_pattern(int index)  {
          return m_results.get_radiation_pattern(index);
  }
          
  /*! \brief Get structure excitation results
  \param index The zero-based index for the result (simulations can return more than one set of results).
  \return The requested radiation pattern data (or NULL if the result does not exist).
  \note You must NOT delete the results object when finished with it.
  */
  inline nec_structure_excitation* get_structure_excitation(int index)  {
          return m_results.get_structure_excitation(index);
  }
  
  /*! \brief Get near field pattern results
  \param index The zero-based index for the result (simulations can return more than one set of results).
  \return The requested radiation pattern data (or NULL if the result does not exist).
  \note You must NOT delete the results object when finished with it.
  */  
  inline nec_near_field_pattern* get_near_field_pattern(int index)  {
          return m_results.get_near_field_pattern(index);
  }
  
  /*! \brief Get structure currents results
  \param index The zero-based index for the result (simulations can return more than one set of results).
  \return The requested radiation pattern data (or NULL if the result does not exist).
  \note You must NOT delete the results object when finished with it.
  */  
  inline nec_structure_currents* get_structure_currents(int index)  {
          return m_results.get_structure_currents(index);
  }
  
  /* added for the python wrapping : some access functions */
  
  void set_isave(int in_isave)
  {
          isave = in_isave;
  }
  
  int get_inc()
  {
          return inc;
  }
  
  nec_float get_xpr1()
  {
          return xpr1;
  }
  
  nec_float get_xpr2()
  {
          return xpr2;
  }
  
  /* end of functions added for the python wrapping */
  
  inline void set_output(nec_output_file in_output, nec_output_flags in_output_flags)
  {
          m_output = in_output;
          m_output_flags = in_output_flags;
          
          m_output_fp = m_output.get_fp();
  }
  
  inline void set_results_format(enum RESULT_FORMAT result_format)
  {
          m_results.m_result_format = result_format;
  }
  
  inline void set_gain_only(bool flag)
  {
          m_output_flags.set_gain_only(flag);
  }
  
  
  /*!\brief Benchmark the libnecpp engine. A score of 100 is roughly an Athlon XP 1800. */
  static nec_float benchmark();
  
  /*! \brief Signal the end of a geometry description.
  
  This function prepares for a calculation by calling calc_prepare().
  */
  void geometry_complete(int gpflag);

  
  /*! Set the prameters of the medium (permittivity and permeability)
  
          \param permittivity The electric permittivity of the medium (in farads per meter)
          \param permeability The magnetic permeability of the medium (in henries per meter)

          From these parameters a speed of light is chosen.
  */
  void medium_parameters(nec_float permittivity, nec_float permeability)  {
    em::constants::permittivity = permittivity;
    em::constants::permeability = permeability;
  }
  
  
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
  void wire(int tag_id, int segment_count,
  nec_float xw1, nec_float yw1, nec_float zw1,
  nec_float xw2, nec_float yw2, nec_float zw2,
  nec_float rad, nec_float rdel, nec_float rrad);

  void sp_card(int ns,
      nec_float x1, nec_float y1, nec_float z1,
      nec_float x2, nec_float y2, nec_float z2);

  void sc_card( int i2,
      nec_float x3, nec_float y3, nec_float z3,
      nec_float x4, nec_float y4, nec_float z4);

  void gx_card(int i1, int i2);

  void move( nec_float rox, nec_float roy, nec_float roz, nec_float xs,
                        nec_float ys, nec_float zs, int its, int nrpt, int itgi );

  /*! Add an arc to the geometry,
  
  All co-ordinates are in meters.
  
    \param tag_id The tag ID.
    \param segment_count The number of segments.
    \param rada The radius.
    \param ang1 The angle of the arc starting point.
    \param ang2 The angle of the arc end point.
    \param rad The wire radius.
  */
  void arc( int tag_id, int segment_count, nec_float rada,
      nec_float ang1, nec_float ang2, nec_float rad );
      
  
  /*! \brief Add an helix to the geometry,
  
  \remark The helix is a versatile geometry element. For example, to generate a spiral printed circuit antenna, use a helix of zero height.

  All co-ordinates are in meters.
  
    \param tag_id The tag ID.
    \param segment_count The number of segments.
    \param s The turn spacing.
    \param h1 The total length of the helix (negative for a left-handed helix).
    \param a1 x-start radius.
    \param b1 y-start radius.
    \param a2 x-end radius.
    \param b2 y-end radius.
    \param rad The wire radius.
  */
  void helix(int tag_id, int segment_count, nec_float s, nec_float hl, nec_float a1, nec_float b1,
      nec_float a2, nec_float b2, nec_float rad);
  
  
  /*! \brief "fr" card, frequency parameters
   \param ifrq \c O= LINEAR STEP \c 1=MULTIPLICATIVE
   \param nfrq Number of steps, for a single frequency use nfrq=1
   \param freq_mhz Frequency or Starting frequency (in MHz)
   \param del_freq Frequency increment, (ADD if ifrq=0, OR MULTIPLY if ifrq=1)
  */
  void fr_card(int ifrq, int nfrq, nec_float freq_mhz, nec_float del_freq);
  
  /*! \brief "ld" card, loading parameters
  \verbatim
  LD  LOADING
    itmp1-   -1 CANCEL LOADS,
      0=SERIES RLC LUMP,
      1=PARALLEL RLC LUMP,
      2=SERIES DIST.,
      3=PARALLEL DIST. (A),
      4=Z (B),
      5=WIRE COND. (C)
    itmp2- TAG# TO BE LOADED, BLANK/0= USE ABSOLUTE #s  
    itmp3- SEG# OF TAG # TO START LOADS, OR ABSOLUTE SEG#
    itmp4- SEG# OF TAG# TO END LOADS, OR OR ABSOLUTE SEG#
    F1- RES., OHMS, OR (A) OHMS/UNIT LENGTH, OR (B) RES. OR (C) OHMS/METER
    F2- IND., HENRY, OR (A) HY/LENGTH OR (B) REACT. OR (C) BLANK
    F3- CAP,. FARAD, OR (A,B) BLANK
  \endverbatim
  */
  void ld_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3);
  
  /*! \brief Ground parameters under the antenna
  
  \remark Specifies the relative dielectric constant and conductivity of ground in the vicinity of the antenna. 
  In addition, a second set of ground parameters for a second medium can be specified, or a radial wire 
  ground screen can be modeled using a reflection coefficient approximation.
  
  \param ground_type (was IPERF) Ground-type flag. The options are:
    \arg \c -1 - Nullifies ground parameters previously used and sets free-space condition. The remainder of the \
    parameters should be zero in this case. 
    \arg \c O - Finite ground, reflection-coefficient approximation.
    \arg \c 1 - Perfectly conducting ground.
    \arg \c 2 - Finite ground, Sommerfeld/Norton method. 
  
  \param rad_wire_count (was NRADL) - Number of radial wires in the ground screen approximation; Set to zero implies no ground screen.
  
  \param EPSE (F1) - Relative dielectric constant for ground in the vicinity of the antenna. Set to zero in case of a perfect ground.
  \param SIG (F2) - Conductivity in mhos/meter of the ground in the vicinity of the antenna. Set to zero in the case of a perfect ground. If SIG is input as a negative number, the complex dielectric constant Ec = Er -j*sigma/(omega*epsilonzero) is set to EPSR - |SIG|.
  
  \remark
   
  Options for Remaining Floating Point Fields (F3-F6):
    \li a. For an infinite ground plane, F3 through F6 are blank. 
    \li b. Radial wire ground screen approximation (NRADL nonzero). The ground screen is always centered at the origin, i.e., at (0,0,0), and lies in the XY plane. (F3) - The radius of the screen in meters. (F4) - Radius of the wires used in the screen, in meters. (F5) & (F6) - Blank.
    \li c. Second medium parameters (NRADL = O) for medium outside the region of the first medium (cliff problem). These parameters alter the far field patterns but do not affect the antenna impedance or current distribution. (F3) - Relative dielectric constant of medium 2. (F4) - Conductivity of medium 2 in mhos/meter. (F5) - Distance in meters from the origin of the coordinate system to the join between medium 1 and 2. This distance is either the radius of the circle where the two media join or the distance out the positive X axis to where the two media join in a line parallel to the Y axis. Specification of the circular or linear option is on the RP card. See Figure 16. (F6) - Distance in meters (positive or zero) by which the surface of medium 2 is below medium 1.
  */
  void gn_card(int ground_type, int rad_wire_count, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
  
  
  /*! "ex" card, excitation parameters
  \verbatim
      EX  EXCITE STRUCTURE, LAST ENCOUNTERED=USED
        I1- 0=E VOLTAGE (A), 1=LINEAR WAVE (B), 2= R CIRC WAVE (B) 3=L CIRC WAVE (B), 4= CURRENT (C), 5= VOLTAGE DISC. (A)
        I2- (A) SOURCE TAG#, (B) # TH ANGLS, (C) BLANK
        I3- (A) SOURCE SEG#, (B) # PH ANGLS, (C) BLANK
        I4- (A) XX= ADMIT.,IMPED. PRINT, X=0 NO/1 DO, (BC), 1= ADM. PRINT
        F1- (A) EREAL, (B) TH ANGL, (C) X OF SOURCE
        F2- (A) EIMAG, (B) PH ANGL, (C) Y OF SOURCE
        F3- (A) NORM FOR I4, (B) ET ANGL, Z OF SOURCE
        F4- (A) BLANK, (B) TH INC, (C) ALPHA ANGLE FROM XY
        F5- (A) BLANK, (B) PH INC, (C) BETA ANGLE FROM X
        F6- (A) BLANK, (B) MIN/MAJ AXIS, PRODUCT AMPS X LENGTH
        
        // NOT YET DONE... F7- (A) BLANK, (B) INCIDENT AMPLITUDE (Volts/m)
  \endverbatim
  */
  void ex_card(enum excitation_type itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);


  /*! \brief Specifies the excitation for the structure. The excitation can be voltage sources on the structure, an elementary current source,
      or a plane-wave incident on the structure.

      All angles are in degrees.

      \param excitation_type Determines the type of excitation which is used :
              excitation_type = O - voltage source (applied-E-field source). 
              excitation_type = 1 - incident plane wave, linear polarization. 
              excitation_type = 2 - incident plane wave, right-hand (thumb along the incident k vector) elliptic polarization. 
              excitation_type = 3 - incident plane wave, left-hand elliptic polarization. 
              excitation_type = 4 - elementary current source. 
              excitation_type = 5 - voltage source (current-slope-discontinuity).

      \param itmp2 If excitation_type = 0 or 5 : the tag number of the source segment  (if itmp1 = 0 absolute segment numbers will be used) ;
                    else if excitation_type = 1, 2 or 3 : number of theta angles desired for the incident plane wave ;
                    else zero.
            
      \param itmp3 If excitation_type = 0 or 5 : the rank (among the segments the tag number of which is itmp2) or absolute segment number
                      of the source segment ;
                    else if excitation_type = 1, 2 or 3 : number of phi angles desired for the incident plane wave ;
                    else zero.
      
      \param itmp4 If itmp4 = 1 the maximum relative admittance matrix asymmetry for source segment (if excitation_type = 0 or 5) and
              network connections (whatever excitation_type may be) will be calculated and printed.
      
      \param itmp5 If excitation_type = 0 or 5 : tmp3 will be taken under account if itmp5 = 1 ;
                    else zero.
      
      \param tmp1 If excitation_type = 0 or 5 : the real part of the voltage  ;
                  else if excitation_type = 1, 2 or 3 : the first value of theta ;
                  else the x-coordinate of the current source.
      
      \param tmp2 If excitation_type = 0 or 5 : the imaginary part of the voltage  ;
                  else if excitation_type = 1, 2 or 3 : the first value of phi ;
                  else if excitation_type = 4 : the y-coordinate of the current source.
      
      \param tmp3 If excitation_type = 0 or 5 : the normalization constant for the impedance printed in the optional impedance table (if tmp3 = 0
                      the impedance will be normalized to their maximum value) ;
                  else if excitation_type = 1, 2 or 3 : eta in degrees. Eta is the polarization angle defined as the angle between the
                      theta unit vector and the direction of the electric field for linear polarization or the major ellipse axis for elliptical polarization ;
                  else if excitation_type = 4 : the z-coordinate of the current source.
      
      \param tmp4 If excitation_type = 0 or 5 : zero.
                  else excitation_type = 1, 2 or 3 : theta angle stepping increment.
                  else if excitation_type = 4 : the angle the current source makes with the XY plane.
          
      \param tmp5 If excitation_type = 0 or 5 : zero.
                  else excitation_type = 1, 2 or 3 : phi angle stepping increment.
                  else if excitation_type = 4 : the angle the projection of the current source on the XY plane makes with the X axis.

      \param tmp6 If excitation_type = 0 or 5 : zero.
                  else excitation_type = 1, 2 or 3 : ratio of minor axis to major axis for elliptic polarization (major axis field strength - 1 V/m).
                  else if excitation_type = 4 : "Current moment" of the source (in amp meter).    
    */
  void ex_card(enum excitation_type itmp1, int itmp2, int itmp3, int itmp4, int itmp5,
                  nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
  {
    int itmp45 = 10*itmp4 + itmp5;
    return this->ex_card( itmp1, itmp2, itmp3, itmp45, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );  
  }
  
  
  
  /*! \brief "tl" card, transmission line parameters
  
  \remark
  To generate a transmission line between any two points on the structure. Characteristic impedance, length, and shunt admittance are the defining parameters.
  
  \verbatim
  TL TRANSMISSION LINE 
    I1- PORT 1 TAG #, BLANK/0, USE I2 AS ABSOLUTE
    I2- SEGMENT#, OR ABSOLUTE END 1 SEGMENT, -1=CANCEL NETS/LINES
    I3- AS I1 FOR PORT 2
    I4- AS I2 FOR PORT 2
    F1- LINE Zo, -=CROSSED LINE
    F2- LINE LENGTH METERS, BLANK=STRAIGHT LINE P1 TO P2
    F3- REAL SHUNT ADM., END 1 MHOS
    F4- IMAG SHUNT ADM., END 1
    F5- REAL SHUNT ADM., END 2
    F6- IMAG SHUNT ADM., END 2
  \endverbatim
  */
  void tl_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
  
  /*! 4: "nt" card, network parameters
  \verbatim
    NT  NETWORKS
      I1- PORT 1 TAG #, BLANK/0, USE I2 AS ABSOLUTE
      I2- SEGMENT#, OR ABSOLUTE END 1 SEGMENT, -1=CANCEL NETS/LINES
      I3- AS I1 FOR PORT 2
      I4- AS I2 FOR PORT 2
      F1- REAL OF Y(11), MHOS
      F2- IMAG OF Y(11)
      F3- REAL OF Y(12)
      F4- IMAG OF Y(12)
      F5- REAL OF Y(22)
      F6- IMAG OF Y(22)
  \endverbatim
  */
  void nt_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
  
  /*! \brief "xq" execute card - calc. including radiated fields
  
  \param itmp1
     \c 0 NO PATTERN,
     \c 1 XY PATTERN,
     \c 2 YZ PATTERN,
     \c 3 BOTH. (DO NOT USE FOR RADIAL GND SCREEN OR 2ND GND MEDIUM)
    
    \remark FOR A SINGLE FREQUENCY, XQ, NE, NH, RP CAUSE IMMEDIATE EXECUTION
    FOR MULTIPLE FREQS, ONLY XQ, RP CAUSE EXECUTION
  */
  void xq_card(int itmp1);
  
  /*! "gd" card, ground representation */
  void gd_card(nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4);
  
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
    
  \li The rp_card() function will call simulate(), causing the interaction matrix to be computed and factored and the structure currents to be computed if these operations have not already been performed. Hence, all required input parameters must be set before the rp_card() function is called. 
  \li At a single frequency, any number of rp_card() calls may occur in sequence so that different field-point spacings may be used over different regions of space. If automatic frequency stepping is being used (i.e., in_nfrq on the fr_card() function is greater than one), only one rp_card() function will act as data inside the loop. Subsequent calls to rp_card() will calculate patterns at the final frequency. 
  \li When both n_theta and n_phi are greater than one, the angle theta (or Z) will be stepped faster than phi. 
  \li When a ground plane has been specified, field points should not be requested below the ground (theta greater than 90 degrees or Z less than zero.)
  
  */
  void rp_card(int calc_mode,
    int n_theta, int n_phi,
    int output_format, int normalization, int D, int A,
    nec_float theta0, nec_float phi0, nec_float delta_theta, nec_float delta_phi,
    nec_float radial_distance, nec_float gain_norm);
  
   /*! "pt" card, print control for current */
  void pt_card(int itmp1, int itmp2, int itmp3, int itmp4);
  
  
   /*! "pq" card, print control for charge */
  void pq_card(int itmp1, int itmp2, int itmp3, int itmp4);
  
  
  
  /*! "kh" card, matrix integration limit */
  void kh_card(nec_float tmp1);
  
  
  /*! Near field calculation parameters 
  
  \remark
  \li If the number of frequencies is not equal to one (as specified by the fr_card() function, then the ne_card() function will call simulate(), causing the interaction matrix to be computed and factored and the structure currents to be computed.
  */
  void ne_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
  
  /*! Near field calculation parameters 
  
  \remark
  \li If the number of frequencies is not equal to one (as specified by the fr_card() function, then the ne_card() function will call simulate(), causing the interaction matrix to be computed and factored and the structure currents to be computed.
  */
  void nh_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
  
  /*! "ek" card,  extended thin wire kernel option */
  void set_extended_thin_wire_kernel(bool ekflag);
  
  
  /*! "cp" card, maximum coupling between antennas */
  void cp_card(int itmp1, int itmp2, int itmp3, int itmp4);
  
  
  /*! "pl" card, plot flags 
    \exception int Throws int on error.
  */
  void pl_card(const char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4);
  
  
  /*!****************************************************
  *** normal exit of nec2++ when all jobs complete ok ***
  ******************************************************/
  void all_jobs_completed()
  {
  }
  
  void write_results(ostream& os)
  {
    m_results.write(os);
  }
  
  
  
  /*!  \brief Start a simulation
  
    This function will trigger a calculation. In the traditional NEC
    world, This signals the end of the main input section and the
    beginning of the frequency do loop.
    
    \param far_field_flag is true if last card was XQ or RP
    \warning far_field_flag is should never be specified as true
    because both the xq_card() and rp_card() functions will call
    this function automatically.
  */
  void simulate(bool far_field_flag = false);

  //! an object to pipe output through...
  nec_output_file m_output;
  nec_ground ground;
  c_geometry* m_geometry;
  c_plot_card plot_card;
  
  c_ground_wave ground_wave;

  
  //! pq card flags
  int iptflq;
  int iptaq, iptaqf, iptaqt;
  
  
  //! pt card flags...
  int iptflg;
  int iptag, iptagf, iptagt;
  
  
  int iflow;
  int ifrq, nfrq;
  nec_float delfrq;
  
  // strcture loading
  int_array ldtyp, ldtag, ldtagf, ldtagt;
  real_array zlr, zli, zlc;
  
  // normalized receiving pattern
  real_matrix fnorm;
  
  int nthi, nphi;
  nec_float thetis, phiss;
  
  
  /*!\brief The results object that holds all the specific results.
  */
  nec_results m_results;
  
  
  //! an object to pipe output through...
  nec_output_flags m_output_flags;
  
  
  nec_float _wavelength;
  
  /* common  /cmb/ */
  complex_array cm;  // primary interaction matrix
  
  /* common  /matpar/ */
  int icase, npblk, nlast;
  int nbbx, npbx, nlbx, nbbl, npbl, nlbl;
  
  /* common  /save/ */
  int_array ip;
  nec_float freq_mhz;
  
  /* common  /crnt/ */
  real_array air, aii; //! coefficients of the constant terms in the current interpolation functions for the current vector 
  real_array bir, bii; //! coefficients of the sine terms in the current interpolation functions
  real_array cir, cii; //! coefficients of the cosine terms in the current interpolation functions
  complex_array current_vector;  //! the current vector
  
  int ifar; //! input integer flag (from RP card) specifies type of field computation, or type of ground system for far fields
    
  /* common  /zload/ */
  int nload;
  complex_array zarray;
  
  /* common  /yparm/ */
  int ncoup, icoup;
  int_array nctag, ncseg;
  complex_array y11a, y12a;
  
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
  
  /* common  /fpat/ */
  enum excitation_type m_excitation_type;
  
  int m_rp_output_format;
  int m_rp_normalization;
  
  int m_near, nfeh, nrx, nry, nrz, nth, nph, ipd, iavp;
  nec_float thets, phis, dth, dph, rfld, gnor; 
  nec_float xpr6, structure_power_loss, xnr, ynr, znr, dxnr, dynr, dznr;
  
  
  /* common  /dataj/ */
  int ind1, indd1, ind2, indd2;
  bool m_use_exk; /* Was iexk */
  
  nec_float m_s, m_b, xj, yj, zj, cabj, sabj, salpj;
  nec_float rkh; /* matrix integration limit */
  nec_float t1xj, t1yj, t1zj, t2xj, t2yj, t2zj;
  nec_complex  exk, eyk, ezk, exs, eys, ezs, exc, eyc, ezc;
  
  /* common  /smat/ */
  int nop; /* My addition */
  complex_array symmetry_array;
  
  /* common  /incom/ */
  int isnor;
  nec_float xo, yo, zo, sn, xsn, ysn;
  
  /* common  /tmi/ */
  int ija; /* changed to ija to avoid conflict */
  nec_float zpk, rkb2;
  
  /*common  /tmh/ */
  nec_float zpka, rhks;

  
  // some auxiliary functions to be made private once
  // the radiation pattern calculation is done entirely
  // inside this class...
  void gfld(nec_float rho, nec_float phi, nec_float rz,
    nec_complex *eth, nec_complex *epi,
    nec_complex *erd, bool space_only, nec_float in_wavelength  );
  
  void ffld(nec_float thet, nec_float phi,
    nec_complex *eth, nec_complex *eph, nec_float in_wavelength  );
  

private:

  /*! \brief A private convenience function called by ne_card() and nh_card()
  */
  void ne_nh_card(int in_nfeh, int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
  
  
  void print_freq_int_krnl(
    nec_float f, 
    nec_float lambda, 
    nec_float int_dist, 
    bool using_extended_kernel);
    
  void   antenna_env(void);
  
  /*no more used */
  void  print_structure_currents(char *pattype, int iptflg, int iptflq, int iptag, int iptagf, int iptagt, int iptaq, 
    int iptaqf, int iptaqt);
    
  /*!\brief Calculate network data such as the lengths of transmission lines.
  */
  void  calculate_network_data(void);
  void  print_network_data(void);
  void  print_norm_rx_pattern();
  void  print_input_impedance();
  void  print_power_budget(void);
  void  structure_segment_loading();
    
  
  enum excitation_return
    excitation_loop(int in_freq_loop_state, int mhz);
      
  void  setup_excitation();
  
  /* pointers to output files */
  FILE *m_output_fp;
  
  int inc, processing_state, isave;
  int nthic, nphic;
  int iped;
  
  nec_float impedance_norm_factor; // was zpnorm
  
  nec_float xpr1, xpr2, xpr3, xpr4, xpr5, xpr7;
  
  nec_structure_currents* structure_currents;
  
  void load();

  void cmset(int64_t nrow, complex_array& in_cm, nec_float rkhx);
  void compute_matrix_ss(int j1, int j2, int im1, int im2,
      complex_array& in_cm, int64_t nrow, int itrp);
  void cmsw(int j1, int j2, int i1, int i2, complex_array& in_cm,
      complex_array& cw, int64_t ncw, int64_t nrow, int itrp);
  void cmws( int j, int i1, int i2, complex_array& in_cm,
                        int64_t nr, complex_array& cw, int64_t nw, int itrp );
        
  void cmww(int j, int i1, int i2, complex_array& in_cm, int64_t nr,
      complex_array& cw, int64_t nw, int itrp);
  void couple(complex_array& cur, nec_float wlam);

  void efld(nec_float xi, nec_float yi, nec_float zi, nec_float ai, bool on_source_segment);
  void eksc(nec_float s, nec_float z, nec_float rh, nec_float xk, int ij,
      nec_complex *ezs, nec_complex *ers, nec_complex *ezc,
      nec_complex *erc, nec_complex *ezk, nec_complex *erk);
  void ekscx(nec_float bx, nec_float s, nec_float z, nec_float rhx, nec_float xk,
      int ij, int inx1, int inx2, nec_complex *ezs,
      nec_complex *ers, nec_complex *ezc, nec_complex *erc,
      nec_complex *ezk, nec_complex *erk);
  void etmns(nec_float p1, nec_float p2, nec_float p3, nec_float p4, nec_float p5,
      nec_float p6, nec_float incident_amplitude, enum excitation_type excite_type, complex_array& e);

  void fblock( int nrow, int ncol, int imax, int ipsym );

  void gf(nec_float zk, nec_float *co, nec_float *si);
  void gh(nec_float zk, nec_float *hr, nec_float *hi);
  void gx(nec_float zz, nec_float rh, nec_float xk,
      nec_complex *gz, nec_complex *gzp);
  void gxx(nec_float zz, nec_float rh, nec_float a, nec_float a2, nec_float xk,
      int ira, nec_complex *g1, nec_complex *g1p, nec_complex *g2,
      nec_complex *g2p, nec_complex *g3, nec_complex *gzp);
  void hfk(nec_float el1, nec_float el2, nec_float rhk,
      nec_float zpkx, nec_float *sgr, nec_float *sgi);
  void hintg(nec_float xi, nec_float yi, nec_float zi);
  void hsfld(nec_float xi, nec_float yi, nec_float zi, nec_float ai);
  void hsflx(nec_float s, nec_float rh, nec_float zpx, nec_complex *hpk,
      nec_complex *hps, nec_complex *hpc);

  void intx(nec_float el1, nec_float el2, nec_float b, int ij,
      nec_float *sgr, nec_float *sgi);

  void nefld(nec_float xob, nec_float yob, nec_float zob, nec_complex *ex,
      nec_complex *ey, nec_complex *ez);
  void netwk(complex_array& in_cm, int_array& in_ip, complex_array& einc);
  void nfpat(void);
  void nhfld(nec_float xob, nec_float yob, nec_float zob, nec_complex *hx,
      nec_complex *hy, nec_complex *hz);
  void pcint(nec_float xi, nec_float yi, nec_float zi, nec_float cabi,
      nec_float sabi, nec_float salpi, complex_array& e);
  void impedance_print(int in1, int in2, int in3, nec_float fl1, nec_float fl2,
      nec_float fl3, nec_float fl4, nec_float fl5, nec_float fl6, const char *ia);
  void qdsrc(int is, nec_complex v, complex_array& e);

  
  void rom2(nec_float a, nec_float b, complex_array& sum, nec_float dmin);
  void sflds(const nec_float t, complex_array& e);
  void solgf(nec_complex *a, nec_complex *b, nec_complex *c,
      nec_complex *d, nec_complex *xy, int *ip, int np, int n1,
      int n, int mp, int m1, int m, int n1c, int n2c, int n2cz);
  void unere(nec_float xob, nec_float yob, nec_float zob, bool ground_reflection);
  nec_complex zint(nec_float sigl, nec_float rolam);
  
  void init_voltage_sources();

}; /* nec_context */

#endif /* __nec_context__ */

