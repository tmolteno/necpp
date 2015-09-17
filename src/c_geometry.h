/*
	Copyright (C) 2004-2011  Timothy C.A. Molteno
	
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
#ifndef __c_geometry__
#define __c_geometry__
#include "math_util.h"
#include <vector>
#include <iostream>

/* Replaces the "10000" limit used to */
/* identify segment/patch connections */
#define	PCHCON  100000

class nec_context;
#include "nec_output.h"
#include "nec_wire.h"



/*! \brief A Class describing the antenna geometry
 * \file c_geometry.h
 * \par Using c_geometry
 * 
  \verbatim
    c_geometry g = new c_geometry();
    g.parse_geometry(file);
 \endverbatim           
  OR
  \verbatim
    c_geometry g = new c_geometry();
    g.wire(xxxxxxx);
    g.arc();
    ...
    g.geometry_complete();
 \endverbatim		
*/
class c_geometry
{
public:
  c_geometry();
  void set_context(nec_context* m_context);
  
  
  /*! \brief Generates segment geometry for a straingt wire
          \param tag_id
          \param segment_count Number of Elements (should be around 12-20 per wavelength)
          \param rad Wire radius of first segment (in Meters)
          \param rdel Ratio of the length of a segment to the length of the previous segment.  (Set to 1.0 if segments have uniform length)
          \param rrad The ratio of the radii of adjacent segments (Set to 1.0 if not tapered)
  */
  void wire( 	int tag_id, int segment_count,
                  nec_float xw1, nec_float yw1, nec_float zw1,	// first co-ordinate
                  nec_float xw2, nec_float yw2, nec_float zw2,	// second co-ordinate
                  nec_float rad,
                  nec_float rdel, nec_float rrad);
                  
                  
  void arc( int tag_id, int segment_count, nec_float rada,
                  nec_float ang1, nec_float ang2, nec_float rad );
                  
  void helix(int tag_id, int segment_count, 
             nec_float s, nec_float hl, nec_float a1, nec_float b1,
             nec_float a2, nec_float b2, nec_float rad);
  
  void move( nec_float rox, nec_float roy, nec_float roz, nec_float xs,
                  nec_float ys, nec_float zs, int its, int nrpt, int itgi );

  /*! \brief Reflects partial structure along x,y, or z axes.
          
    \param ix If ix = 1 then the structure is reflected along X axis. 
    \param iy If iy = 1 then the structure is reflected along Y axis.
    \param iz If iz = 1 then the structure is reflected along Z axis.
    \param itx The tag number increment.    
  */
  void reflect(int ix, int iy, int iz, int itx)  {
    int nop = 100*ix + 10*iy + iz;
    this->reflect(ix, iy, iz, itx, nop);
  }

  /*! \brief Rotates structure along Z-axis to generate a copies in a cylindrical array.

    \param itx The tag number increment.
    \param nop The total number of times that the structure is to occur in the cylindrical array.  
  */
  void generate_cylindrical_structure(int itx, int nop)
  {
    this->reflect(-1, 0, 0, itx, nop);
  }

  void reflect( int ix, int iy, int iz, int itx, int nop );

  /*! \brief Scale all dimensions of a structure by a constant.*/
  void scale( nec_float xw1);
  
  void patch( int nx, int ny,
                nec_float ax1, nec_float ay1, nec_float az1,
                nec_float ax2, nec_float ay2, nec_float az2,
                nec_float ax3, nec_float ay3, nec_float az3,
                nec_float ax4, nec_float ay4, nec_float az4 );
  
  
  void sp_card( int ns,
                nec_float in_x1, nec_float in_y1, nec_float in_z1,
                nec_float in_x2, nec_float in_y2, nec_float in_z2);
  
  void sc_card( int i2,
                nec_float x3, nec_float y3, nec_float z3,
                nec_float x4, nec_float y4, nec_float z4);
  void sc_multiple_card(int i2,
        nec_float x3, nec_float y3, nec_float z3,
        nec_float x4, nec_float y4, nec_float z4);
  
  void gx_card(int card_int_1, int card_int_2);
       
  /*! \brief Geometry is complete
   *  \exception nec_exception* If there is an error with the geometry.
   */
  void geometry_complete(nec_context* m_context, int gpflag);
  
  
  
  /*!\brief Parse an NEC geometry description contained in the file input_fp
  */
  void parse_geometry(nec_context* m_context, FILE* input_fp);

  
  /*!\brief Helper method to decide whether extended. thin-wire approximation can be used 
  */
  int test_ek_approximation(int seg1, int seg2);
  
  int get_segment_number( int in_tag, int m);

  
  void frequency_scale(nec_float freq_mhz);

  void tbf( int i, int icap);
  void trio( int j );
  
  void get_current_coefficients(nec_float wavelength, complex_array& curx,
                  real_array& air, real_array& aii,
                  real_array& bir, real_array& bii,
                  real_array& cir, real_array& cii,
                  complex_array& vqds, int nqds,
                  int_array& iqds);
  
  nec_float patch_angle(int patch_index, nec_float in_ax, nec_float in_ay, nec_float in_az);
  
  /*! \brief Calculate the xyz components of the electric field due to surface currents.
  */
  void fflds(nec_float rox, nec_float roy, nec_float roz,
          complex_array& scur, 
          nec_complex *in_ex, nec_complex *in_ey, nec_complex *in_ez );
          
  int n_segments;	// The number of segments
  int np;
  int_array segment_tags;
  real_array x, y, z, segment_length, segment_radius;
  real_array x2, y2, z2;
  real_array cab, sab, salp;
  
  int m, mp;	// The number of patches
  int m_ipsym;
  
  real_array t1x, t1y, t1z, t2x, t2y, t2z;	// t1, t2 basis co-ordinates?
  real_array px, py, pz, pbi, psalp;			// patch data

  
  int_array icon1, icon2;
  // Connected Segment Information (was common  /segj/  in FORTRAN code)
  int jsno, nscon, maxcon; /* Max. no. connections */
  int_array jco;
  real_array ax, bx, cx;
  
  inline int n_plus_m(void) const {
    return n_segments + m;
  }
  int n_plus_2m, n_plus_3m; /* n+m,n+2m,n+3m */
private:
  //	The geometry data measured in meters is stored in these arrays
  //	and the x,y,z,si,bi arrays are then scaled for each frequency
  real_array x_unscaled, y_unscaled, z_unscaled, si_unscaled, bi_unscaled;
  real_array px_unscaled, py_unscaled, pz_unscaled, pbi_unscaled;
  
  void sbf( int i, int is, nec_float *aa, nec_float *bb, nec_float *cc );

  void divide_patch( int nx );

  void connect_segments( int ignd );

  void read_geometry_card(FILE* input_fp, char *gm,
          int *i1, int *i2, 
          nec_float *x1, nec_float *y1,nec_float *z1,
          nec_float *x2, nec_float *y2, nec_float *z2, 
          nec_float *rad );

  nec_context* m_context;
  nec_output_file* m_output;

  std::vector<nec_wire> m_wires;
  
  int patch_type;
  nec_3vector patch_x1, patch_x2, patch_x3, patch_x4;
  bool   _prev_sc;
};


#endif /* __c_geometry__ */
