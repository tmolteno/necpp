/*
  Copyright (C) 2004-2015  Timothy C.A. Molteno
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
#include "nec_context.h"
#include "c_geometry.h"
#include "nec_exception.h"

#include <cstdlib>

nec_context::nec_context() : fnorm(0,0), current_vector(0) {
  m_output_fp=NULL;
  m_geometry = new c_geometry();
  inc=0;
  isave=0;
  nthic=0;
  nphic=0;
  
  impedance_norm_factor = 0.0; // was zpnorm
  xpr1=0.0;
  xpr2=0.0;
  xpr3=0.0;
  xpr4=0.0;
  xpr5=0.0;
  xpr7=0.0;
  
  /*structure_currents is a pointer to the "nec_base_result" which takes care of storing and printing of the currents*/
  structure_currents = NULL;
          
  // allocate the ground grid
  initialize();
}

nec_context::~nec_context() {
  delete m_geometry;
}

/*! \brief
Initialize everything, called after construction so that we
can tell the geometry object what nec_context object to
point to.
*/
void nec_context::initialize() {
  DEBUG_TRACE("initialize()");
  nthi=0;
  nphi=0;
  iflow = 1;
  thetis=0.0;
  phiss=0.0;
  
  iptag=0;
  iptagf=0;
  iptagt=0;
  iptaq=0;
  iptaqf=0;
  iptaqt=0;
  
  init_voltage_sources();
  m_geometry->set_context(this);
}

/*! \brief Private method to initialize the voltage source buffers
*/
void nec_context::init_voltage_sources() {
  /* Free vsource buffers */
  ivqd.resize(0);
  iqds.resize(0);
  vqd.resize(0);
  vqds.resize(0);
  source_segment_array.resize(0);
  source_voltage_array.resize(0);

  voltage_source_count=0;
  nvqd=0;
  iped=0;
}

/*! \brief After the geometry has been specified, this function prepares for calculations
*/
void nec_context::calc_prepare()
{
  DEBUG_TRACE("calc_prepare()");
  iflow=1;    
  
  int n_plus_m = m_geometry->n_plus_m();
  /* Allocate some buffers */
  air.resize(n_plus_m);
  aii.resize(n_plus_m);
  bir.resize(n_plus_m);
  bii.resize(n_plus_m);
  cir.resize(n_plus_m);
  cii.resize(n_plus_m);
  
  ip.resize(m_geometry->n_plus_2m);
        ASSERT((m_geometry->n_segments + m_geometry->m*3 == m_geometry->n_plus_3m));
        // TODO FIX HERE... THIS SIZE AINT RIGHT
  current_vector.resize(m_geometry->n_segments + m_geometry->m*4);
  
  /* Matrix parameters */
  neq= m_geometry->n_plus_2m;
  neq2=0;

  /* default values for input parameters and flags */
  npeq = m_geometry->np + 2*m_geometry->mp;
  processing_state=1;
  rkh=1.;
  m_use_exk=false;
  m_excitation_type = EXCITATION_VOLTAGE;
  nload=0;
  network_count=0;
  m_near = -1;
  ifar = -1;
  ncoup=0;
  icoup=0;
  freq_mhz= em::speed_of_light() / 1.0e6;
  ground.default_values();
  nfrq=1;
  iptflg = -2;
  iptflq = -1;
  iped=0;
}


/*!\brief Benchmark the libnecpp engine. A score of 1.0 is roughly an Athlon XP 1800. */
nec_float nec_context::benchmark()
{
  nec_float start_timer, stop_timer;
  
  secnds( &start_timer );
  for (int i=0; i<20; i++)
  {
    {
      /*
      CMEXAMPLE 2. CENTER FED LINEAR ANTENNA. 
      CM           CURRENT SLOPE DISCONTINUITY SOURCE. 
      CM           1. THIN PERFECTLY CONDUCTING WIRE 
      CE           2. THIN ALUMINUM WIRE 
      GW 0 8 0. 0. -.25 0. 0. .25 .00001 
      GE 
      FR 0 3 0 0 200. 50. 
      EX 5 0 5 1 1. 0. 50. 
      XQ 
      LD 5 0 0 0 3.720E+07 
      FR 0 1 0 0 300. 
      EX 5 0 5 0 1. 
      GN 1
      XQ 
      EN
      */

      nec_context nec;
      nec.set_gain_only(true);
      nec.initialize();
      
      c_geometry* geo = nec.get_geometry();
      geo->wire(0,8, 0.0, 0.0, -0.25,
          0.0, 0.0, 0.25,
          0.00001,
          1.0, 1.0);
          
      nec.geometry_complete(0);
      
      nec.fr_card(0, 3, 200.0, 50.0);
      nec.ex_card(EXCITATION_VOLTAGE_DISC, 0, 5, 1, 1.0, 0.0, 50.0, 0.0, 0.0, 0.0);
      nec.xq_card(0);
      nec.ld_card(5, 0, 0,0, 3.72e7, 0.0, 0.0);
      nec.fr_card(0, 1, 300.0, 0.0);
      nec.ex_card(EXCITATION_VOLTAGE_DISC, 0, 5, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      nec.gn_card(1, 0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0);
      
      // TODO Keep testing here...
      nec.xq_card(0);
  /*  This line kills the version on Wine! Very strange! I think that the gn card is bad
      fixme:msvcrt:MSVCRT_signal (11 (nil)):stub
      err:seh:EXC_DefaultHandling Unhandled exception code c0000005 flags 0 addr 0x404e2380
      Wine failed with return code 1
  */
    }
    {
      /*
        An example showing how to evaluate for maximum gain.
        
      CM GA - NEC FILE
      CE go blue ! 
      GW 0 36 0 0 0 -0.042 0.008 0.017 0.001
      GW 0 21 -0.042 0.008 0.017 -0.048 0.021 -0.005 0.001
      GW 0 70 -0.048 0.021 -0.005 0.039 0.032 -0.017 0.001
      GW 0 70 -0.048 0.021 -0.005 0.035 0.043 0.014 0.001
      GW 0 50 -0.042 0.008 0.017 0.017 -0.015 0.014 0.001
      GW 0 66 0.017 -0.015 0.014 -0.027 0.04 -0.031 0.001
      GW 0 85 -0.027 0.04 -0.031 0.046 -0.01 0.028 0.001
      GW 0 47 0.046 -0.01 0.028 -0.013 -0.005 0.031 0.001
      GW 0 70 0.017 -0.015 0.014 -0.048 -0.038 -0.04 0.001
      GW 0 77 -0.048 -0.038 -0.04 0.049 -0.045 -0.04 0.001
      GE 0
      GN -1
      LD 5 0 0 0 3.720E+07 
      FR 0 1 0 0 2400
      PT -1
      EX 1 1 1 0 0 0 0 0 0 0 0
      RP 0 1 1 0500 90 90 0 0
      EN
      
      */
      nec_context nec;
      nec.set_gain_only(false);
      nec.initialize();
      
      c_geometry* geo = nec.get_geometry();
      geo->wire(0, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
      geo->wire(0, 21, 0.042, 0.002, 0.017, -0.028, 0.011, 0.005, 0.001, 1.0, 1.0);
      geo->wire(0, 70, -0.058, 0.021, 0.005, 0.039, 0.062, 0.017, 0.001, 1.0, 1.0);
      geo->wire(0, 70, 0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0);
      geo->wire(0, 50, 0.042, 0.018, 0.017, 0.017, 0.015, 0.024, 0.001, 1.0, 1.0);
      geo->wire(0, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0);
      geo->wire(0, 85, 0.027, 0.04, 0.031, -0.046, -.01, 0.028, 0.001, 1.0, 1.0);
      geo->wire(0, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0);
      geo->wire(0, 70, 0.027, 0.015, 0.014, -0.078, 0.038, -0.04, 0.001, 1.0, 1.0);
      geo->wire(0, 77, 0.078, -0.038, 0.04, 0.049, 0.065, 0.04, 0.001, 1.0, 1.0);
      nec.geometry_complete(0);
      
      nec.gn_card(-1,0,0.0, 0.0, 0.0,0.0, 0.0, 0.0);
      nec.ld_card(5,0,0,0,3.72e7,0.0,0.0);
      nec.pt_card(-1, 0, 0, 0);
      nec.ex_card(EXCITATION_LINEAR, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      nec.fr_card(0, 2, 2400.0, 100.0);
      nec.rp_card(0, 1, 1, 0,5,0,0, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0);
      nec.get_gain_max();
    }
  }
  /* time the process */
  secnds( &stop_timer );
  stop_timer -= start_timer;
  cout << endl << endl;
  nec_float sec = stop_timer / 1000.0;
  
  nec_float bench = 70.0 / sec;
  
  return bench;
}


/*! \brief Signal the end of a geometry description.

This function prepares for a calculation by calling calc_prepare().
*/
void nec_context::geometry_complete(int gpflag) {
  DEBUG_TRACE("geometry_complete()");
  m_geometry->geometry_complete(this, gpflag);
  calc_prepare();
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
void nec_context::wire(int tag_id, int segment_count,
  nec_float xw1, nec_float yw1, nec_float zw1,
  nec_float xw2, nec_float yw2, nec_float zw2,
  nec_float rad, nec_float rdel, nec_float rrad)
{
  m_geometry->wire(tag_id, segment_count, xw1, yw1, zw1, xw2, yw2, zw2, rad, rdel, rrad);
}


void nec_context::sp_card(int ns,
    nec_float x1, nec_float y1, nec_float z1,
    nec_float x2, nec_float y2, nec_float z2)
{
  m_geometry->sp_card(ns, x1, y1, z1, x2, y2, z2);
}



void nec_context::sc_card(int i2,
    nec_float x3, nec_float y3, nec_float z3,
    nec_float x4, nec_float y4, nec_float z4)
{
  m_geometry->sc_card(i2, x3, y3, z3, x4, y4, z4);
}

void nec_context::gx_card(int i1, int i2) {
  m_geometry->gx_card(i1, i2);
}

void nec_context::move( nec_float rox, nec_float roy, nec_float roz, nec_float xs,
    nec_float ys, nec_float zs, int its, int nrpt, int itgi )
{
  DEBUG_TRACE("nec_context::move nrpt=" << nrpt);
  m_geometry->move(rox, roy, roz, xs, ys, zs, its, nrpt, itgi);
}


/*! Add an arc to the geometry,

All co-ordinates are in meters.

  \param tag_id The tag ID.
  \param segment_count The number of segments.
  \param rada The radius.
  \param ang1 The angle of the arc starting point.
  \param ang2 The angle of the arc end point.
  \param rad The wire radius.
*/
void nec_context::arc( int tag_id, int segment_count, nec_float rada,
        nec_float ang1, nec_float ang2, nec_float rad ) {
    m_geometry->arc(tag_id, segment_count, rada, ang1, ang2, rad);
}
    

/*! \brief Add an helix to the geometry,

\remark The helix is a versatile m_geometry->element. For example, to generate a spiral printed circuit antenna, use a helix of zero height.

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
void nec_context::helix(int tag_id, int segment_count, nec_float s, nec_float hl, nec_float a1, nec_float b1,
    nec_float a2, nec_float b2, nec_float rad) {
    
  /* int tag_id, int segment_count, 
             nec_float s, nec_float hl, nec_float a1, nec_float b1,
             nec_float a2, nec_float b2, nec_float rad
  */
  m_geometry->helix(tag_id, segment_count, s, hl, a1, b1, a2, b2, rad);
}

/* "fr" card, frequency parameters

FREQUENCY
I1- O= LINEAR STEP, 1=MULTIPLICATIVE
I2- NO. STEPS, BLANK=1
I3- BLANK -- not used in this function
I4- BLANK -- not used in this function
F1- FREQUENCY OR START FREQUENCY
F2- FREQ INCREMENT, ADD OR MULTIPLY
*/
void nec_context::fr_card(int in_ifrq, int in_nfrq, nec_float in_freq_mhz, nec_float in_del_freq) {
  DEBUG_TRACE("fr_card()");
  ifrq = in_ifrq;
  nfrq = in_nfrq;
  if ( nfrq == 0)
          nfrq=1;
          
  freq_mhz = in_freq_mhz;
  delfrq = in_del_freq;
  if ( iped == 1)
          impedance_norm_factor = 0.0;
          
  processing_state = 1;
  iflow = 1;
}



void nec_context::ld_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3)
{
  DEBUG_TRACE("ld_card()");
  if ( iflow != 3 )
  {
    iflow=3;
    /* Free loading buffers */
    nload=0;
    ldtyp.resize(0);
    ldtag.resize(0);
    ldtagf.resize(0);
    ldtagt.resize(0);
    zlr.resize(0);
    zli.resize(0);
    zlc.resize(0);
  
    if ( processing_state > 2 )
      processing_state=2;
    if ( itmp1 == -1 )
      return; // continue card input loop
  }

  /* Reallocate loading buffers */
  nload++;
  ldtyp.resize(nload);
  ldtag.resize(nload);
  ldtagf.resize(nload);
  ldtagt.resize(nload);
  zlr.resize(nload);
  zli.resize(nload);
  zlc.resize(nload);

  int idx = nload-1;
  ldtyp[idx]= itmp1;
  ldtag[idx]= itmp2;
  if ( itmp4 == 0)
    itmp4= itmp3;
  ldtagf[idx]= itmp3;
  ldtagt[idx]= itmp4; 

  if ( itmp4 < itmp3 )
  {
    nec_stop("DATA FAULT ON LOADING CARD No: %d: ITAG "
      "STEP1: %d IS GREATER THAN ITAG STEP2: %d",
      nload, itmp3, itmp4 );
  }

  zlr[idx]= tmp1;
  zli[idx]= tmp2;
  zlc[idx]= tmp3;
}



/* "gn" card, ground parameters under the antenna
GN    NEAR GROUND, GROUND SCREEN, ADDED GROUND
  I1- -1=SET FREE SPACE (A), 0=REFL COEFF, 1=IDEAL (B), 2-SOMMERFELD
  I2- (A) BLANK), NO WIRES IN GND SCREEN (C), 0= NO WIRES (D)
  I3- BLANK
  I4- BLANK
  F1- (A,B) BLANK, DIELECTRIC OF NEAR GROUND
  F2- (A,B) BLANK, CONDUCTIVITY OF NEAR GROUND
  F3- (A,B) BLANK, (C) RADIUS OF SCREEN, (D) DIELECTRIC 2ND MEDIUM
  F4- (A,B) BLANK, (C) RADII SCREEN WIRES, (D) CONDUCT. 2ND MEDIUM
  F5- (A,B) BLANK, (C) BLANK, (D) DIST TO 2ND MEDIUM, SEE RP
  F6- (A,B) BLANK, (C) BLANK, (D) HEIGHT 2ND MEDIUM (AS IN GD)
*/
void nec_context::gn_card(int ground_type, int rad_wire_count, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  DEBUG_TRACE("gn_card(" << ground_type << ")");
  ground.parse_gn(ground_type, rad_wire_count, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
  
  iflow=4;

  if ( processing_state > 2)
    processing_state=2;
}



/*
    EX  EXCITE STRUCTURE, LAST ENCOUNTERED=USED
      I1- 0=E VOLTAGE (A), 1=LINEAR WAVE (B), 2= R CIRC WAVE (B)
      3=L CIRC WAVE (B), 4= CURRENT (C), 5= VOLTAGE DISC. (A)
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
*/
void nec_context::ex_card(enum excitation_type itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  DEBUG_TRACE("ex_card(" << itmp1
    << "," << itmp2 
    << "," << itmp3 
    << "," << itmp4 
    << "," << tmp1 
    << "," << tmp2 
    << "," << tmp3 
    << "," << tmp4 
    << "," << tmp5 
    << "," << tmp6
    << ")");

  if ( iflow != 5)  {
    init_voltage_sources();

    iflow=5;
    
    if ( processing_state > 3)
      processing_state=3;
  }

  masym = itmp4/10;
  ASSERT(itmp1 >= 0);
  
  m_excitation_type = itmp1;
  if ( (m_excitation_type == EXCITATION_VOLTAGE) || (m_excitation_type == EXCITATION_VOLTAGE_DISC) )  // (Voltage excitation)
  {
    ntsol=0;
  
    if ( m_excitation_type == EXCITATION_VOLTAGE_DISC) // Voltage DISC. 
    {
      nvqd++;
      ivqd.resize(nvqd );
      iqds.resize(nvqd );
      vqd.resize( nvqd );
      vqds.resize(nvqd );
    
      int indx = nvqd-1;
  
      ivqd[indx]= m_geometry->get_segment_number( itmp2, itmp3);
      vqd[indx]= nec_complex( tmp1, tmp2);
      if ( abs( vqd[indx]) < 1.e-20)
        vqd[indx] = cplx_10();
  
      iped= itmp4- masym*10;
      impedance_norm_factor= tmp3;
      if ( (iped == 1) && (impedance_norm_factor > 0.0) )
        iped=2;
      return; /* continue card input loop */
    } /* if ( m_excitation_type == EXCITATION_VOLTAGE_DISC) */
  
    voltage_source_count++;
    source_segment_array.resize(voltage_source_count );
    source_voltage_array.resize(voltage_source_count );
  
    {
      int indx = voltage_source_count-1;
    
      int seg_number = m_geometry->get_segment_number( itmp2, itmp3);
      if (seg_number > m_geometry->segment_length.size())  {
        nec_exception* nex = new nec_exception("CHECK DATA, PARAMETER SPECIFYING EXCITATION SOURCE SEGMENT [");
        nex->append(seg_number);
        nex->append("] IS TOO LARGE" );
        throw nex;
      }
      source_segment_array[indx] = seg_number;
      
      DEBUG_TRACE("Voltage Source: " << nec_complex( tmp1, tmp2));
      source_voltage_array[indx]= nec_complex( tmp1, tmp2);
      if ( abs( source_voltage_array[indx]) < 1.e-20)
        source_voltage_array[indx] = cplx_10();
    
      iped= itmp4- masym*10;
      impedance_norm_factor= tmp3;
      if ( (iped == 1) && (impedance_norm_factor > 0.0) )
        iped = 2;
      return; /* continue card input loop */
    }  
  } /* if ( (m_excitation_type == 0) || (m_excitation_type == 5) ) */

  nthi= itmp2;
  nphi= itmp3;
  xpr1= tmp1;
  xpr2= tmp2;
  xpr3= tmp3;
  xpr4= tmp4;
  xpr5= tmp5;
  xpr6= tmp6;
  // xpr7= tmp7; Put this in here once we are parsing NEC4 excitation stuff.
  voltage_source_count=0;
  nvqd=0;
  thetis= xpr1;
  phiss= xpr2;
}


/* 5: "tl" cards, network parameters

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

*/
void nec_context::tl_card(int itmp1, int itmp2, int itmp3, int itmp4,
  nec_float characteristic_impedance, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  if ( iflow != 6)  {
    network_count=0;
    ntsol=0;
    iflow=6;
  
    if ( processing_state > 3)
      processing_state=3;
  
    if ( itmp2 == -1 )
      return; /* continue card input loop */
  }

  /* Re-allocate network buffers */
  network_count++;
  ntyp.resize(network_count);
  iseg1.resize(network_count);
  iseg2.resize(network_count);
  x11r.resize(network_count);
  x11i.resize(network_count);
  x12r.resize(network_count);
  x12i.resize(network_count);
  x22r.resize(network_count);
  x22i.resize(network_count);

  int idx = network_count-1;
  ntyp[idx] = 2; // TL card

  iseg1[idx]= m_geometry->get_segment_number( itmp1, itmp2);
  iseg2[idx]= m_geometry->get_segment_number( itmp3, itmp4);
  x11r[idx]= characteristic_impedance;
  x11i[idx]= tmp2;
  x12r[idx]= tmp3;
  x12i[idx]= tmp4;
  x22r[idx]= tmp5;
  x22i[idx]= tmp6;

  if (characteristic_impedance <= 0.0)  {
    // Negative characteristic impedance implies a crossed line
    ntyp[idx] = 3;
    x11r[idx] = -characteristic_impedance;
  }
    
}

/* 4:
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

*/
void nec_context::nt_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  if ( iflow != 6)  {
    network_count=0;
    ntsol=0;
    iflow=6;
  
    if ( processing_state > 3)
      processing_state=3;
  
    if ( itmp2 == -1 )
      return; /* continue card input loop */
  }

  /* Re-allocate network buffers */
  network_count++;
  ntyp.resize(network_count);
  iseg1.resize(network_count);
  iseg2.resize(network_count);
  x11r.resize(network_count);
  x11i.resize(network_count);
  x12r.resize(network_count);
  x12i.resize(network_count);
  x22r.resize(network_count);
  x22i.resize(network_count);

  int idx = network_count-1;
  ntyp[idx]=1; // NT card

  iseg1[idx]= m_geometry->get_segment_number( itmp1, itmp2);
  iseg2[idx]= m_geometry->get_segment_number( itmp3, itmp4);
  x11r[idx]= tmp1;
  x11i[idx]= tmp2;
  x12r[idx]= tmp3;
  x12i[idx]= tmp4;
  x22r[idx]= tmp5;
  x22i[idx]= tmp6;
}

/* "xq" execute card - calc. including radiated fields

XQ  EXECUTE ACCUMULATED CARD DECK
  itmp1-
    0=NO PATTERN,
    1=XY PATTERN,
    2= YZ PATTERN,
    3=BOTH
  (DO NOT USE FOR RADIAL GND SCREEN OR 2ND GND MEDIUM)

  NOTES: FOR A SINGLE FREQUENCY, XQ, NE, NH, RP CAUSE IMMEDIATE EXECUTION
  FOR MULTIPLE FREQS, ONLY XQ, RP CAUSE EXECUTION
*/
void nec_context::xq_card(int itmp1)
{
  DEBUG_TRACE("xq_card(" << itmp1 << ")");
  DEBUG_TRACE("iflow =" << iflow);
  if (   ((iflow == 10) && (itmp1 == 0)) ||
    ((nfrq  ==  1) && (itmp1 == 0) && (iflow > 7)) )
    return; /* continue card input loop */

  if ( itmp1 == 0)  {
    if ( iflow > 7)
      iflow=11;
    else
      iflow=7;
  } else {
    ifar=0;
    rfld=0.;
    ipd=0;
    iavp=0;
    m_rp_normalization=0;
    m_rp_output_format=0;
    nth=91;
    nph=1;
    thets=0.0;
    phis=0.0;
    dth=1.0;
    dph=0.0;
  
    if ( itmp1 == 2)
      phis=90.0;
  
    if ( itmp1 == 3)
    {
      nph=2;
      dph=90.0;
    }
  } /* if ( itmp1 == 0) */
  
  simulate(true);
}

/* "gd" card, ground representation */
void nec_context::gd_card(nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4) {
  DEBUG_TRACE("gd_card(" << tmp1 << ")");
  ground.setup_cliff(tmp1, tmp2, tmp3, tmp4);
  iflow=9;
}

/*! \brief Standard radiation pattern parameters 

\param calc_mode
\param n_theta
\param n_phi
\param output_format The output format (0 major axis, minor axis and total gain printed. 
    1 vertical, horizontal ant total gain printed.)
\param normalization Controls the type of normalization of the radiation pattern
N = 0 no normalized gain. 
= 1 major axis gain normalized. 
= 2 minor axis gain normalized. 
= 3 vertical axis gain normalized. 
= 4 horizontal axis gain normalized. 
= 5 total gain normalized.

\param D Selects either power gain or directive gain for both standard printing and normalization. If the structure excitation is an incident plane wave, the quantities printed under the heading "gain" will actually be the scattering cross section (a/lambda 2 ) and will not be affected by the value of d. The column heading for the output will still read "power" or "directive gain," however. 
<ul>
<li> D = 0 power gain. 
<li> D = 1 directive gain.
</ul>

\param A - Requests calculation of average power gain over the region covered by field points. 
<ul>
<li>A = 0 no averaging. 
<li>A = 1 average gain computed. 
<li>A = 2 average gain computed, printing of gain at the field points used for averaging is suppressed. If NTH or NPH is equal to one, average gain will not be computed for any value of A since the area of the region covered by field points vanishes.
</ul>

\param theta0 - Initial theta angle in degrees (initial z coordinate in meters if I1 = 1).

\param phi0 - Initial phi angle in degrees.

\param delta_theta - Increment for theta in degrees (increment for z in meters if I1 = 1).

\param delta_phi - Increment for phi in degrees.

\param radial_distance - Radial distance (R) of field point from the origin in meters. radial_distance is optional. If it is zero, the radiated electric field will have the factor exp(-jkR)/R omitted. If a value of R is specified, it should represent a point in the far-field region since near components of the field cannot be obtained with an RP card. (If I1 = 1, then radial_distance represents the cylindrical coordinate phi in meters and is not optional. It must be greater than about one wavelength.)

\param gain_norm - Determines the gain normalization factor if normalization has been requested in the normalization parameter. If gain_norm is zero, the gain will be normalized to its maximum value. If gain_norm is not zero, the gain w111 be normalized to the value of gain_norm.
*/
void nec_context::rp_card(int calc_mode, int n_theta, int n_phi,
  int output_format, int normalization, int D, int A,
  nec_float theta0, nec_float phi0,
  nec_float delta_theta, nec_float delta_phi,
  nec_float radial_distance, nec_float gain_norm)
{
  DEBUG_TRACE("rp_card(" << calc_mode << ")");
  ifar= calc_mode;
  nth = n_theta;
  nph = n_phi;

  if ( nth == 0)
    nth=1;
  if ( nph == 0)
    nph=1;

  m_rp_output_format = output_format;
  m_rp_normalization = normalization;
  ipd = D;
  iavp = A;
  
  DEBUG_TRACE(" xnda = (" << m_rp_output_format << m_rp_normalization << ipd << iavp << ")");
  if ( m_rp_output_format != 0)
    m_rp_output_format=1;
  if ( ipd != 0)
    ipd=1;
    
  // sanity check. Not point normalizing if there are too few data points to normalize on
  if ( (nth < 2) || (nph < 2) || (ifar == 1) )
    iavp=0;

  thets = theta0;
  phis = phi0;
  dth = delta_theta;
  dph = delta_phi;
  rfld = radial_distance;
  gnor = gain_norm;
  iflow=10;
  
  simulate(true);
}

  /* "pt" card, print control for current */
void nec_context::pt_card(int itmp1, int itmp2, int itmp3, int itmp4) {
  iptflg= itmp1;
  iptag= itmp2;
  iptagf= itmp3;
  iptagt= itmp4;

  if ( (itmp3 == 0) && (iptflg != -1) )
    iptflg = -2;
  if ( itmp4 == 0)
    iptagt= iptagf;
    
}


  /* "pq" card, print control for charge */
void nec_context::pq_card(int itmp1, int itmp2, int itmp3, int itmp4) {
  iptflq= itmp1;
  iptaq= itmp2;
  iptaqf= itmp3;
  iptaqt= itmp4;

  if ( (itmp3 == 0) && (iptflq != -1) )
    iptflq = -2;
  if ( itmp4 == 0)
    iptaqt= iptaqf;
}



/* "kh" card, matrix integration limit */
void nec_context::kh_card(nec_float tmp1) {
  rkh = tmp1;
  if ( processing_state > 2)
    processing_state=2;
  iflow=1;
}

void nec_context::ne_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  ne_nh_card(0, itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

void nec_context::nh_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  ne_nh_card(1, itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

/* \brief Near field calculation parameters 
*/
void nec_context::ne_nh_card(int in_nfeh, int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6)
{
  nfeh = in_nfeh;

  if ( (iflow == 8) && (nfrq != 1) ) {
    m_output.endl(2);
    m_output.line("WHEN MULTIPLE FREQUENCIES ARE REQUESTED, "
            "ONLY ONE NEAR FIELD CARD CAN BE USED -");
    m_output.line("  LAST CARD READ WILL BE USED" );
  }

  m_near= itmp1;
  nrx= itmp2;
  nry= itmp3;
  nrz= itmp4;
  xnr= tmp1;
  ynr= tmp2;
  znr= tmp3;
  dxnr= tmp4;
  dynr= tmp5;
  dznr= tmp6;
  iflow=8;

  if ( nfrq == 1)
    simulate();
}


/* "ek" card,  extended thin wire kernel option */
void nec_context::set_extended_thin_wire_kernel(bool ek_flag) {
  m_use_exk = ek_flag;
  
  if ( processing_state > 2)
    processing_state=2;
  iflow=1;  
}


/* "cp" card, maximum coupling between antennas */
void nec_context::cp_card(int itmp1, int itmp2, int itmp3, int itmp4) {
  if ( iflow != 2) {
    ncoup=0;
    nctag.resize(0);
    ncseg.resize(0);
    y11a.resize(0);
    y12a.resize(0);
  }

  icoup=0;
  iflow=2;

  if ( itmp2 == 0)
    return; /* continue card input loop */

  ncoup++;
  nctag.resize(ncoup);
  ncseg.resize(ncoup);
  nctag[ncoup-1]= itmp1;
  ncseg[ncoup-1]= itmp2;

  if ( itmp4 == 0)
    return; /* continue card input loop */

  ncoup++;
  nctag.resize(ncoup);
  ncseg.resize(ncoup);
  nctag[ncoup-1]= itmp3;
  ncseg[ncoup-1]= itmp4;
}


/* "pl" card, plot flags 
  throws int on error.
*/
void nec_context::pl_card(const char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4) {
  std::string fname(ploutput_filename);
  plot_card = c_plot_card(itmp1,itmp2,itmp3,itmp4, fname);
}



/*! \brief Start a simulation

    This function will trigger a calculation. In the traditional NEC
    world, This signals the end of the main input section and the
    beginning of the frequency do loop.
    
    \param far_field_flag is true if last card was XQ or RP
    \warning far_field_flag is should never be specified as true
    because both the xq_card() and rp_card() functions will call
    this function automatically.
*/
void nec_context::simulate(bool far_field_flag) {
  DEBUG_TRACE("simulate(" << far_field_flag << ")");
  
  /* Allocate the normalization buffer */
  if ( iped )
    fnorm.resize(nfrq,4);
  if ( iptflg >= 2 )
    fnorm.resize(nthi, nphi);


  /* igox is a state variable that is used to change from
    one processing mode to another. The processing to be 
    performed are as follows:
      
    1: Memory allocation & Initialization
    2: Structure segment loading
    3: Excitation set up (right hand side, -e inc.)
    4: ?
    5: Near field calculation
    6: standard far field calculation
  */
  int igox;
  int mhz = 0;
  
  if ( (far_field_flag == true)
    && (processing_state == 5) )
    igox = 6;
  else
    igox = processing_state;
  
  
  try  {
  int64_t iresrv = 0;
  bool in_freq_loop = false;
  
  do {  
    switch( igox )
    {
    case 1: /* Memory allocation for primary interacton matrix. */

      if (false == in_freq_loop)  {
        // TODO fix up the following (changed 2* to 3*)
        iresrv = (m_geometry->n_plus_2m) * (m_geometry->np+3*m_geometry->mp);
        cm.resize(iresrv);
      
        /* Memory allocation for symmetry array */
        nop = neq/npeq;
        symmetry_array.resize(nop*nop);  
        mhz = 1;
        
        fblock( npeq, neq, iresrv, m_geometry->m_ipsym);
        
        in_freq_loop = true;
      }
      
      if ( mhz != 1) {
        if ( ifrq == 1)
          freq_mhz *= delfrq;
        else
          freq_mhz += delfrq;
      }

      _wavelength = em::get_wavelength(1.0e6 * freq_mhz);

      print_freq_int_krnl(freq_mhz, _wavelength, rkh, m_use_exk);
        
      m_geometry->frequency_scale(freq_mhz);
      processing_state = 2;
      /* Falls through. */

    case 2: /* structure segment loading */
      structure_segment_loading();

      processing_state=3;
      ntsol=0;
      /* Falls through. */
    
    case 3: /* excitation set up (right hand side, -e inc.) */
      nthic=1;
      nphic=1;
      inc=1;
      nprint=0;
      /* Falls through. */

    default:
      enum excitation_return ret = excitation_loop(igox, mhz);
    
      if (FREQ_LOOP_CONTINUE == ret)
      {
        continue; // Continue frequency loop
      }
      if (FREQ_LOOP_CARD_CONTINUE == ret)
      {
        throw 1; // Continue card input
      }
    
      nphic = 1;
  
      /* normalized receiving pattern printed */
      print_norm_rx_pattern();
      
      xpr2  = phiss;
  
      if ( mhz == nfrq)
        ifar = -1;
  
      if ( nfrq == 1)
      {
        m_output.end_section();
        throw 1; // Continue card input
      }
    
      print_input_impedance();
    
      nfrq=1;
      mhz=1;    
    } /* switch( igox ) */
  }
  while( (++mhz <= nfrq) );
  } /* try */
  catch (int excep)
  {
    ASSERT(excep == 1);
    // keep going on the card input. The exception
    // is thrown in order to continue card input.
  }
}

/* ********************************************************************************************************** */




void nec_context::print_freq_int_krnl(
  nec_float f, 
  nec_float lambda, 
  nec_float int_dist, 
  bool using_extended_kernel)
{
  m_output.end_section();
  m_output.set_indent(31);
  m_output.line("--------- FREQUENCY --------");
  m_output.string("FREQUENCY= "); m_output.real_out(11,4,f); m_output.line(" MHZ");
  m_output.string("WAVELENGTH="); m_output.real_out(11,4,lambda); m_output.line(" METERS");

  m_output.endl(2);
  m_output.set_indent(24);
  m_output.line("APPROXIMATE INTEGRATION EMPLOYED FOR SEGMENTS");
  m_output.string("THAT ARE MORE THAN "); m_output.real_out(5,3,int_dist,false); m_output.line(" WAVELENGTHS APART");

  if ( using_extended_kernel )
    m_output.line( "THE EXTENDED THIN WIRE KERNEL WILL BE USED");
    
  m_output.set_indent(0);
}

void nec_context::antenna_env(void) {
  ground.calculate_antenna_environment(ground_wave, freq_mhz);
  ground.output_antenna_environment(m_output);
}


/*!\brief Calculate network data such as the lengths of transmission lines.
  TODO Fix up this horrible mess. As far as I can figure we should simply
  be calculating the x11i[] value for every network with ntyp = 2 or 3 (transmission lines).
  The bug being reported by Remi, could be caused by this network failing to calculate the
  length correctly under some circumstances.

      // This is an extremely strange statement.
      //  ntyp[j]    net_type  ntyp[j]/net_type
      //  1      1        1
      //  2      1        2
      //  3      1        3
      //  1      2        0
      //  2      2        1
      //  3      2        1
*/
void nec_context::calculate_network_data(void)
{
  if ( (network_count == 0) || (inc > 1) )
    return;
    
  int itmp3 = 0;
  int net_type = ntyp[0];

  for (int i = 0; i < 2; i++ ) {
    if ( net_type == 3)
      net_type = 2;

    for (int j = 0; j < network_count; j++) {
      if ( (ntyp[j]/net_type) != 1 ) {
        itmp3 = ntyp[j]; // can never be zero
      } else {
        if ( (ntyp[j] >= 2) && (x11i[j] <= 0.0) ) {
          int idx4 = iseg1[j]-1;
          int idx5 = iseg2[j]-1;
          nec_float xx = m_geometry->x[idx5]- m_geometry->x[idx4];
          nec_float yy = m_geometry->y[idx5]- m_geometry->y[idx4];
          nec_float zz = m_geometry->z[idx5]- m_geometry->z[idx4];
          
          // set the length of the transmission line to be the 
          // straight line distance.
          x11i[j] = _wavelength*sqrt(xx*xx + yy*yy + zz*zz);
        }
      }
    }
    if ( itmp3 == 0)  // can only be zero if all the networks are the same type
      return;
    net_type = itmp3;
  }
} /* calculate_network_data */


void nec_context::print_network_data(void)
{
  int itmp1, itmp2, itmp3, itmp4, itmp5;
  const char *pnet[3] = { "        ", "STRAIGHT", " CROSSED" };
  
  if ( (network_count != 0) && (inc <= 1) )
  {
    m_output.nec_printf( "\n\n\n"
        "                                            "
        "---------- NETWORK DATA ----------" );

    itmp3=0;
    itmp1= ntyp[0];

    for(int i = 0; i < 2; i++ )
    {
      if ( itmp1 == 3)
        itmp1=2;

      if ( itmp1 == 2)
      {
        m_output.endl();
        m_output.nec_printf(
            "  -- FROM -  --- TO --      "
            "TRANSMISSION LINE       "
            " --------- SHUNT ADMITTANCES (MHOS) "
            "---------   LINE\n"
            "  TAG   SEG  TAG   SEG    IMPEDANCE      "
            "LENGTH    "
            " ----- END ONE -----      "
            "----- END TWO -----   TYPE\n"
            "  No:   No:  No:   No:         OHMS      "
            "METERS      REAL      IMAGINARY      "
            "REAL      IMAGINARY" );
      } else if (itmp1 == 1) {
        m_output.endl();
        m_output.nec_printf(
            "  -- FROM -  --- TO --            "
            "--------"
            " ADMITTANCE MATRIX ELEMENTS (MHOS) "
            "---------\n"
            "  TAG   SEG  TAG   SEG   "
            "----- (ONE,ONE) ------  "
            " ----- (ONE,TWO) -----   "
            "----- (TWO,TWO) -------\n"
            "  No:   No:  No:   No:      REAL      "
            "IMAGINARY     "
            " REAL     IMAGINARY       REAL      "
            "IMAGINARY" );
      }
      for (int j = 0; j < network_count; j++) {
        itmp2= ntyp[j];

        if ( (itmp2/itmp1) != 1 ) {
          itmp3 = itmp2;
        } else {
          int idx4, idx5;

          itmp4= iseg1[j];
          itmp5= iseg2[j];
          idx4 = itmp4-1;
          idx5 = itmp5-1;
          m_output.endl();
          m_output.nec_printf(
            " %4d %5d %4d %5d  "
            "%11.4E %11.4E  %11.4E %11.4E  "
            "%11.4E %11.4E %s",
            m_geometry->segment_tags[idx4], itmp4, m_geometry->segment_tags[idx5], itmp5,
              x11r[j], x11i[j], x12r[j], x12i[j],
              x22r[j], x22i[j], pnet[itmp2-1] );
        } /* if (( itmp2/ itmp1) == 1) */
      } /* for( j = 0; j < network_count; j++) */
      if ( itmp3 == 0)
        break;
      itmp1= itmp3;
    } /* for( j = 0; j < network_count; j++) */
  } /* if ( (network_count != 0) && (inc <= 1) ) */
} /* print_network_data */

void nec_context::print_norm_rx_pattern()
{
  if ((iptflg != 2)&&(iptflg != 3))
    return; // do not print
    
  {
    // Call the new nec_results class that handles
    // a normalized receiving pattern.
    
    nec_float theta_step = xpr4;
    nec_float phi_step = xpr5;
    nec_float eta = xpr3;
    nec_float axial_ratio = xpr6;
    string pol_type(nec_structure_currents::hpol(m_excitation_type));
    int segment_number = isave;
    
    nec_norm_rx_pattern* rx_pattern = new nec_norm_rx_pattern(nthi, nphi,
      fnorm,
      thetis, theta_step,
      phiss, phi_step,
      eta, 
      axial_ratio, 
      segment_number, 
      pol_type);
      
    rx_pattern->set_frequency(freq_mhz/(1.e-6));
    m_results.add(rx_pattern);
    
    // write the restuls to our output.
    std::stringstream ss;
    rx_pattern->write_to_file(ss);
    m_output.line(ss.str().c_str());
  }
  
#if 0
  nec_float phi = phiss;
  int itmp1 = nthi * nphi;

  nec_float norm_factor = fnorm[0];
  
  for(int j = 1; j < itmp1; j++ )
    if ( fnorm[j] > norm_factor)
      norm_factor = fnorm[j];

  m_output.endl(3);
  m_output.nec_printf(
    "                     "
    "---- NORMALIZED RECEIVING PATTERN ----\n"
    "                      "
    "NORMALIZATION FACTOR: %11.4E\n"
    "                      "
    "ETA: %7.2f DEGREES\n"
    "                      "
    "TYPE: %s\n"
    "                      AXIAL RATIO: %6.3f\n"
    "                      SEGMENT No: %d\n\n"
    "                      "
    "THETA     PHI       ---- PATTERN ----\n"
    "                      "
    "(DEG)    (DEG)       DB     MAGNITUDE",
    norm_factor, xpr3, hpol[m_excitation_type-1], xpr6, isave );
  
  for(int j = 0; j < nphi; j++ )
  {
    int itmp2 = nthi*j;

    for(int i = 0; i < nthi; i++ )
    {
      int itmp3 = i + itmp2;

      if ( itmp3 < itmp1)
      {
        nec_float _tmp2 = fnorm[itmp3] / norm_factor;
        nec_float _tmp3 = db20( _tmp2);

        m_output.endl();
        m_output.nec_printf("                    %7.2f  %7.2f   %7.2f  %11.4E",
          xpr1, phi, _tmp3, _tmp2 );

        xpr1 += xpr4;
      }

    } /* for( i = 0; i < nthi; i++ ) */

    xpr1= thetis;
    phi += xpr5;
  } /* for( j = 0; j < nphi; j++ ) */
#endif
}

void nec_context::print_power_budget(void)
{
  if ( (m_excitation_type == EXCITATION_VOLTAGE) || (m_excitation_type == EXCITATION_VOLTAGE_DISC) )
  {
    nec_float  radiated_power = input_power- network_power_loss - structure_power_loss;
    nec_float  efficiency = 100.0 * radiated_power/input_power;

    m_output.endl(3);
    m_output.nec_printf(
        "                               "
        "---------- POWER BUDGET ---------\n"
        "                               "
        "INPUT POWER   = %11.4E Watts\n"
        "                               "
        "RADIATED POWER= %11.4E Watts\n"
        "                               "
        "STRUCTURE LOSS= %11.4E Watts\n"
        "                               "
        "NETWORK LOSS  = %11.4E Watts\n"
        "                               "
        "EFFICIENCY    = %7.2f Percent",
        input_power, radiated_power, structure_power_loss, network_power_loss, efficiency );
  } /* if ( (m_excitation_type == EXCITATION_VOLTAGE) || (m_excitation_type == EXCITATION_VOLTAGE_DISC) ) */
}

void nec_context::print_input_impedance() {
  nec_float tmp1, tmp2, tmp3, tmp4, tmp5;

  if ( iped != 0) {
    int iss;

    if ( nvqd > 0)
      iss = ivqd[nvqd-1];
    else
      iss = source_segment_array[voltage_source_count-1];

    m_output.endl(3);
    m_output.nec_printf(
        "                            "
        " -------- INPUT IMPEDANCE DATA --------\n"
        "                                     "
        " SOURCE SEGMENT No: %d\n"
        "                                  "
        " NORMALIZATION FACTOR:%12.5E\n\n"
        "              ----------- UNNORMALIZED IMPEDANCE ----------  "
        "  ------------ NORMALIZED IMPEDANCE -----------\n"
        "      FREQ    RESISTANCE    REACTANCE    MAGNITUDE    PHASE  "
        "  RESISTANCE    REACTANCE    MAGNITUDE    PHASE\n"
        "       MHz       OHMS         OHMS         OHMS     DEGREES  "
        "     OHMS         OHMS         OHMS     DEGREES",
        iss, impedance_norm_factor );

    if ( 0 == ifrq )
      tmp1= freq_mhz-( nfrq-1)* delfrq;
    else
      tmp1= freq_mhz/( pow(delfrq, (nfrq-1)) );

    for(int i = 0; i < nfrq; i++ ) {
      tmp2= fnorm(i,0)/ impedance_norm_factor;
      tmp3= fnorm(i,1)/ impedance_norm_factor;
      tmp4= fnorm(i,2)/ impedance_norm_factor;
      tmp5= fnorm(i,3);

      m_output.endl();
      m_output.nec_printf(
          " %9.3f   %11.4E  %11.4E  %11.4E  %7.2f  "
          " %11.4E  %11.4E  %11.4E  %7.2f",
          tmp1, fnorm(i,0), fnorm(i,1), fnorm(i,2),
          fnorm(i,3), tmp2, tmp3, tmp4, tmp5 );

      if ( ifrq == 0)
        tmp1 += delfrq;
      else
        tmp1 *= delfrq;
    } /* for( i = 0; i < itmp1; i++ ) */
    m_output.end_section();
  } /* if ( iped != 0) */
}


void nec_context::structure_segment_loading()
{
  nec_float tim1, tim, tim2;

  m_output.end_section();
  m_output.line("                          ------ STRUCTURE IMPEDANCE LOADING ------" );

  if ( nload != 0)
    load();

  if ( nload == 0 )
  {  
    m_output.line("                                 THIS STRUCTURE IS NOT LOADED" );
  }
  antenna_env();

  /* fill and factor primary interaction matrix */
  secnds( &tim1 );
  cmset( neq, cm, rkh );
  secnds( &tim2 );
  tim= tim2- tim1;
  factrs(m_output, npeq, neq, cm, ip );
  secnds( &tim1 );
  tim2= tim1- tim2;
  
  m_output.end_section();
  m_output.line("                             ---------- MATRIX TIMING ----------");
  m_output.string("                               FILL= ");
  m_output.integer(int(tim)); 
  m_output.string(" msec  FACTOR: "); 
  m_output.integer(int(tim2)); m_output.string(" msec");
}


void nec_context::setup_excitation()
{
  nec_float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

  tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=0.0;

  if ( (m_excitation_type != 0) && (m_excitation_type != 5) )
  {
    if ( (iptflg <= 0) || (m_excitation_type == 4) )
    {
      m_output.endl(3);
      m_output.line("                             ---------- EXCITATION ----------" );
    }
    
    tmp5= degrees_to_rad(xpr5);
    tmp4= degrees_to_rad(xpr4);

    if ( m_excitation_type == 4)
    {
      tmp1= xpr1/ _wavelength;
      tmp2= xpr2/ _wavelength;
      tmp3= xpr3/ _wavelength;
      tmp6= xpr6/( _wavelength* _wavelength);

      m_output.endl();
      m_output.line(    "                                      CURRENT SOURCE");
      m_output.line(    "                     -- POSITION (METERS) --       ORIENTATION (DEG)");
      m_output.line(    "                     X          Y          Z       ALPHA        BETA   DIPOLE MOMENT");
      m_output.nec_printf("               %10.5f %10.5f %10.5f  %7.2f     %7.2f    %8.3f",
        xpr1, xpr2, xpr3, xpr4, xpr5, xpr6 );
    }
    else
    {
      tmp1= degrees_to_rad(xpr1);
      tmp2= degrees_to_rad(xpr2);
      tmp3= degrees_to_rad(xpr3);
      tmp6= xpr6;

      if ( iptflg <= 0)
      {
        m_output.endl();
        m_output.nec_printf(
          "PLANE WAVE - THETA: %7.2f deg, PHI: %7.2f deg,"
          " ETA=%7.2f DEG, TYPE - %s  AXIAL RATIO: %6.3f",
          xpr1, xpr2, xpr3, nec_structure_currents::hpol(m_excitation_type).c_str(), xpr6 );
      }
    } /* if ( m_excitation_type == 4) */
  }

  nec_float incident_amplitude = 1.0; // the incident electric field amplitude...
  if (xpr7 != 0.0)
    incident_amplitude = xpr7;
    

  /* Fill E-field right-hand matrix */
  etmns( tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, incident_amplitude, m_excitation_type, current_vector);
}

enum excitation_return nec_context::excitation_loop(int in_freq_loop_state, int mhz)
{
  int itmp1;
  
  do
  {
    if (in_freq_loop_state < 4)
    {
      setup_excitation();

      /* matrix solving  (netwk calls solves) */
      calculate_network_data();
      print_network_data();

      if ( (inc > 1) && (iptflg > 0) )
        nprint=1;

      // c_network::net_solve( cm, &cm[ib11], &cm[ic11], &cm[id11], ip, current_vector );
      netwk( cm, ip, current_vector );
      ntsol=1;

      if ( iped != 0)
      {
        itmp1= ( mhz-1);

        fnorm(itmp1,0) = real( zped);
        fnorm(itmp1,1) = imag( zped);
        fnorm(itmp1,2) = abs( zped);
        fnorm(itmp1,3) = arg_degrees( zped);

        if ( iped != 2 )
        {
            if ( fnorm(itmp1,2) > impedance_norm_factor)
             impedance_norm_factor = fnorm(itmp1,2);
        }
      } /* if ( iped != 0) */

      /* printing structure currents */
      
      /* The printing of currents is now managed by a dedicated nec_base_result : nec_structure_currents*/
      
      /* Check that the flag for currents and charge densities don't have unexpected values */ 
      if((iptflg <-2)||(iptflg > 3))
      {
        m_output.line("Warning : The print control flag for currents was uncorrect ; it has been set to -1 (no printing for currents).");
        iptflg = -1;
      }
      
      if((iptflq <-2)||(iptflq > 0))
      {
        m_output.line("Warning : The print control flag for charge densities was uncorrect ; it has been set to -1 (no printing for charge densities).");
        iptflq = -1;
      }
      
      /* Check for compatibility between the excitation type and the output format for currents */      
      if(m_excitation_type != 1 && m_excitation_type != 2 && m_excitation_type != 3)
      {
        if (iptflg == 2 || iptflg == 3)
        {
          m_output.line("Warning : The output format chosen for currents was incompatible with the excitation type. No currents has been printed.");
          iptflg = -1;
        }
      }
      
      /* Now proceed with the printing of currents */      
      if (iptflg != -1 || iptflq != -1)
      {
        if((iptflg < 1 || iptflg == 3 || ((iptflg == 1 || iptflg == 2) && inc <= 1)) || (iptflg == -1 && iptflq != -1))
        {
          structure_currents = new nec_structure_currents(this, m_excitation_type,
          nload, xpr3, xpr6);
        
          if(iptflg != 3)
          {
            structure_currents->set_frequency(freq_mhz/(1.e-6));
            m_results.add(structure_currents);
          }
        }
        structure_currents->analyze();
        
        if(iptflg != 3)
        {  
          std::stringstream ss;
          structure_currents->write_to_file(ss);
          m_output.string(ss.str().c_str(), false);
        }
      
      } /* if (iptflg != -1 || iptflq != -1) */
                
      
      /*print_structure_currents(hpol[m_excitation_type-1], iptflg, iptflq,
        iptag, iptagf, iptagt, iptaq, iptaqf, iptaqt);*/
            
      print_power_budget();

      processing_state = 4;

      if ( ncoup > 0)
        couple( current_vector, _wavelength );

      if ( iflow == 7)
      {
        if ( (m_excitation_type > 0) && (m_excitation_type < 4) )
        {
          nthic++;
          inc++;
          xpr1 += xpr4;

          if ( nthic <= nthi )
            continue; /* continue excitation loop */

          nthic=1;
          xpr1= thetis;
          xpr2= xpr2+ xpr5;
          nphic++;

          if ( nphic <= nphi )
            continue; /* continue excitation loop */

          return FREQ_PRINT_NORMALIZATION;
        } /* if ( (m_excitation_type >= 1) && (m_excitation_type <= 3) ) */

        if ( nfrq != 1) 
        {
          return FREQ_LOOP_CONTINUE; /* continue the freq loop */
        }

        m_output.end_section();
        return FREQ_LOOP_CARD_CONTINUE; /* continue card input loop */

      } /*if ( iflow == 7) */

    }
    
    if (in_freq_loop_state < 5)
      processing_state = 5;

    
    /* near field calculation */
    if (in_freq_loop_state < 6)
    {
      if ( m_near != -1)
      {
        if((nrx*nry*nrz) != 0)
          nfpat();

        if ( mhz == nfrq)
          m_near = -1;

        if ( nfrq == 1)
        {
          m_output.end_section();
            return FREQ_LOOP_CARD_CONTINUE; /* continue card input loop */
        }
      } /* if ( near != -1) */
    }
    
    /* Standard far-field calculation */
    
    if ( ifar != -1)
    {
      nec_radiation_pattern* rad_pat =
        new nec_radiation_pattern(nth, nph,
            thets, phis, dth, dph,
            rfld, ground,
            ifar, _wavelength,
            input_power, network_power_loss,
            m_rp_output_format, m_rp_normalization, ipd, iavp,
            gnor, plot_card);

      rad_pat->analyze(this);
      
      /* Here we write the largest gain to the standard output
        if the user has specified this on the command line.
      */
      if (m_output_flags.get_gain_flag())
      {
        rad_pat->write_gain_normalization();
        delete rad_pat;
      }
      else
      {
        rad_pat->set_frequency(freq_mhz/(1.e-6));
        m_results.add(rad_pat);
        
        // write the results to the ordinary NEC output file.
        
        std::stringstream ss;
        rad_pat->write_to_file(ss);
        m_output.line(ss.str().c_str());
        // print_radiation_pattern(input_power, network_power_loss);
      }
    }

    if ( (m_excitation_type == 0) || (m_excitation_type >= 4) )
    {
      if ( mhz == nfrq )
        ifar = -1;

      if ( nfrq != 1)
      {
        return FREQ_LOOP_CONTINUE;  /* continue the freq loop */
      }

      m_output.end_section();
      return FREQ_LOOP_CARD_CONTINUE;  /* continue card input loop */
    }

    nthic++;
    inc++;
    xpr1 += xpr4;

    if ( nthic <= nthi )
      continue; /* continue excitation loop */

    nthic = 1;
    xpr1  = thetis;
    xpr2 += xpr5;
    nphic++;

    if ( nphic > nphi )
      return FREQ_PRINT_NORMALIZATION;

  }
  while( true );
} /* excitation_loop */



/* load calculates the impedance of specified */
/* segments for various types of loading */
void nec_context::load()
{
  int istep, istepx, l1, l2, ldtags, ichk;
  bool iwarn = false;
  nec_complex zt, tpcj;
  
  int n = m_geometry->n_segments;
  int np = m_geometry->np;
  
  tpcj = nec_complex(0.0,1.883698955e+9);
  m_output.endl();
  m_output.line("  LOCATION        RESISTANCE  INDUCTANCE  CAPACITANCE     IMPEDANCE (OHMS)   CONDUCTIVITY  CIRCUIT");
  m_output.line("  ITAG FROM THRU     OHMS       HENRYS      FARADS       REAL     IMAGINARY   MHOS/METER      TYPE");
  
  /* Initialize impedance array, used for temporary */
  /* storage of loading information. */
  zarray.resize(n);
  zarray.setConstant(cplx_00());
  
  istep = 0;
  
  /* Surely This should be rewritten as...
    for (int istepx=0; istepx<nload; istepx++)
    {
      if ( ldtyp[istepx] > 5 )
      {
        m_output.nec_printf(
          "\n  IMPROPER LOAD TYPE CHOSEN,"
          " REQUESTED TYPE IS %d", ldtyp[istepx] );
        stop(-1);
      }
      ... rest of code here
    }
    // done with loads.
    if ( iwarn == true )
      m_output.nec_printf(
        "\n  NOTE, SOME OF THE ABOVE SEGMENTS "
        "HAVE BEEN LOADED TWICE - IMPEDANCES ADDED" );
  
    if ( nop != 1)
    {
      for( i = 0; i < np; i++ )
      {
        zt= zarray[i];
        l1= i;
      
        for( l2 = 1; l2 < nop; l2++ )
        {
          l1 += np;
          zarray[l1]= zt;
        }
      }
    }
    
    
  */
  /* cycle over loading cards */
  while( true )
  {
    istepx = istep;
    istep++;
  
    if ( istep > nload)
    {
      if ( iwarn == true )
        m_output.line("NOTE, SOME OF THE ABOVE SEGMENTS HAVE BEEN LOADED TWICE - IMPEDANCES ADDED" );
    
      if ( nop == 1)
        return;
    
      for(int i = 0; i < np; i++ )
      {
        zt= zarray[i];
        l1= i;
      
        for( l2 = 1; l2 < nop; l2++ )
        {
          l1 += np;
          zarray[l1]= zt;
        }
      }
      return;
    
    } /* if ( istep > nload) */
  
    if ( ldtyp[istepx] > 5 )
    {
      nec_stop("IMPROPER LOAD TYPE CHOSEN,"
        " REQUESTED TYPE IS %d", ldtyp[istepx] );
    }
  
    /* search segments for proper itags */
    ldtags = ldtag[istepx];
    int jump = ldtyp[istepx]+1;
    ichk=0;
    l1= 1;
    l2= n;
  
    if ( ldtags == 0)
    {
      if ( (ldtagf[istepx] != 0) || (ldtagt[istepx] != 0) )
      {
        l1= ldtagf[istepx];
        l2= ldtagt[istepx];
      }
    }
  
    for(int i = l1-1; i < l2; i++ )
    {
      if ( ldtags != 0)
      {
        if ( ldtags != m_geometry->segment_tags[i])
          continue;
      
        if ( ldtagf[istepx] != 0)
        {
          ichk++;
          if ( (ichk < ldtagf[istepx]) || (ichk > ldtagt[istepx]) )
            continue;
        }
        else
          ichk=1;
      } /* if ( ldtags != 0) */
      else
        ichk=1;
    
      nec_float seg_length = m_geometry->segment_length[i];

      ASSERT(0.0 != m_geometry->segment_length[i]);
      ASSERT(0.0 != m_geometry->segment_radius[i]);
      
      /*   Calculation of lamda*imped. per unit length,
        jump to appropriate section for loading type */
      switch( jump )
      {
        case 1:
          zt= zlr[istepx]/ seg_length+ tpcj* zli[istepx]/( seg_length*_wavelength);
          if ( fabs( zlc[istepx]) > 1.0e-20)
            zt += _wavelength/( tpcj* seg_length* zlc[istepx]);
          break;
      
        case 2:
          zt= tpcj* seg_length* zlc[istepx]/ _wavelength;
          if ( fabs( zli[istepx]) > 1.0e-20)
            zt += seg_length* _wavelength/( tpcj* zli[istepx]);
          if ( fabs( zlr[istepx]) > 1.0e-20)
            zt += seg_length/ zlr[istepx];
          zt=1./ zt;
          break;
      
        case 3:
          zt= zlr[istepx]* _wavelength+ tpcj* zli[istepx];
          if ( fabs( zlc[istepx]) > 1.0e-20)
            zt += 1./( tpcj* seg_length* seg_length* zlc[istepx]);
          break;
      
        case 4:
          zt= tpcj* seg_length* seg_length* zlc[istepx];
          if ( fabs( zli[istepx]) > 1.0e-20)
            zt += 1./( tpcj* zli[istepx]);
          if ( fabs( zlr[istepx]) > 1.0e-20)
            zt += 1./( zlr[istepx]* _wavelength);
          zt=1./ zt;
          break;
      
        case 5:
          zt= nec_complex( zlr[istepx], zli[istepx])/ seg_length;
          break;
      
        case 6:
        {
          zt= zint( zlr[istepx]* _wavelength, m_geometry->segment_radius[i]);
        }
      
      } /* switch( jump ) */
    
      if (( fabs( real( zarray[i]))+ fabs( imag( zarray[i]))) > 1.0e-20)
        iwarn=true;
      zarray[i] += zt;
    
    } /* for( i = l1-1; i < l2; i++ ) */
    
    if ( ichk == 0 )
    {
      nec_exception* nex = new nec_exception("LOADING DATA CARD ERROR, NO SEGMENT HAS AN ITAG = ");
      nex->append(ldtags);
      throw nex;
    }
  
    /* printing the segment loading data, jump to proper print */
    switch( jump )
    {
      case 1:
        impedance_print( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
          zli[istepx], zlc[istepx],0.,0.,0.," SERIES ");
        break;
    
      case 2:
        impedance_print( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
          zli[istepx], zlc[istepx],0.,0.,0.,"PARALLEL");
        break;
    
      case 3:
        impedance_print( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
          zli[istepx], zlc[istepx],0.,0.,0., "SERIES (PER METER)");
        break;
    
      case 4:
        impedance_print( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
          zli[istepx], zlc[istepx],0.,0.,0.,"PARALLEL (PER METER)");
        break;
    
      case 5:
        impedance_print( ldtags, ldtagf[istepx], ldtagt[istepx],0.,0.,0.,
          zlr[istepx], zli[istepx],0.,"FIXED IMPEDANCE ");
        break;
    
      case 6:
        impedance_print( ldtags, ldtagf[istepx], ldtagt[istepx],
          0.,0.,0.,0.,0., zlr[istepx],"  WIRE  ");
    } /* switch( jump ) */
  
  } /* while( true ) */
}

/* cmset sets up the complex structure matrix in the array in_cm */
void nec_context::cmset( int64_t nrow, complex_array& in_cm, nec_float rkhx) {
  int mp2, it, i1, i2, in2;
  int im1, im2, ist, ij, jss, jm1, jm2, jst;
  complex_array scm;
  
  int np = m_geometry->np;
  int mp = m_geometry->mp;
  
  mp2 = 2 * mp;
//   = np + mp2;
//  neq= m_geometry->n_plus_2m; // n + 2* m;
  
  rkh= rkhx;
  //iout=2* npblk* nrow;
  it= nlast;
                  
  vector_fill(in_cm,0,it*nrow,cplx_00());
  
  i1= 1;
  i2= it;
  
  in2= i2; 
  if ( in2 > np)
    in2= np;
  
  
  im1= i1- np;
  im2= i2- np;
  
  if ( im1 < 1)
    im1=1;
  
  ist=1;
  if ( i1 <= np)
    ist= np- i1+2;

  ASSERT(in2 == std::min(nlast, np));
  
  /* wire source loop */
  for( int j = 1; j <= m_geometry->n_segments; j++ ) {
    m_geometry->trio(j);
    
    for (int i = 0; i < m_geometry->jsno; i++ ) {
      ij = m_geometry->jco[i];
      m_geometry->jco[i] = ((ij-1)/ np)* mp2+ij;
    }

    if ( i1 <= in2)
      cmww( j, i1, in2, in_cm, nrow, in_cm, nrow,1);

    if ( im1 <= im2) {
      complex_array temp = in_cm.segment((ist-1)*nrow, in_cm.size() - ((ist-1)*nrow));
      cmws( j, im1, im2, temp, nrow, in_cm, nrow, 1);
      /* CMWS (J,IM1,IM2,CM(1,IST),NROW,CM,NROW,1) */
    }
    /* matrix elements modified by loading */
    if ( nload == 0)
      continue;

    if ( j > np)
      continue;

    int ipr = j;
    if ( (ipr < 1) || (ipr > it) )
      continue;

    nec_complex zaj= zarray[j-1];

    for (int i = 0; i < m_geometry->jsno; i++ ) {
      jss = m_geometry->jco[i];
      in_cm[(jss-1)+(ipr-1)*nrow] -= ( m_geometry->ax[i]+ m_geometry->cx[i])* zaj;
    }
  } /* for( j = 1; j <= n; j++ ) */
  
  int m = m_geometry->m;
  if ( m != 0)  {
    /* matrix elements for patch current sources */
    jm1=1- mp;
    jm2=0;
    jst=1- mp2;

    for (int i = 0; i < nop; i++ ) {
      jm1 += mp;
      jm2 += mp;
      jst += npeq;

      if ( i1 <= in2) {
        int tmp_n = (jst-1);
        complex_array temp = in_cm.segment(tmp_n, in_cm.size()-tmp_n);
        cmsw( jm1, jm2, i1, in2, temp, in_cm, 0, nrow, 1);
      }

      if ( im1 <= im2) {
        int tmp_n = (jst-1)+(ist-1)*nrow;
        complex_array temp = in_cm.segment(tmp_n, in_cm.size()-tmp_n);
        compute_matrix_ss( jm1, jm2, im1, im2, temp, nrow, 1);
      }
    }

  } /* if ( m != 0) */
  
  if ( icase == 1)
    return;
  
  /* Allocate to scratch memory */
  scm.resize(m_geometry->n_plus_2m);
  
  /* combine elements for symmetry modes */
  for (int i = 0; i < it; i++ ) {
    int row_offset = i*nrow;
    
    for( int j = 0; j < npeq; j++ ) {
      for( int k = 0; k < nop; k++ ) {
        int ka = j+ k*npeq;
        scm[k] = in_cm[row_offset + ka];
      }

      in_cm[row_offset + j] = scm.segment(0,nop).sum(); //scm.sum(0,nop);
      
      for( int k = 1; k < nop; k++ ) {
        int ka = j+ k*npeq;
        nec_complex deter = scm[0];

        for( int kk = 1; kk < nop; kk++ ) {
          deter += scm[kk] * symmetry_array[k+kk*nop];
        }
        in_cm[row_offset + ka] = deter;
      }
    }
  }
  
  scm.resize(0);
}

/* compute_matrix_ss computes matrix elements for surface-surface interactions. */

void nec_context::compute_matrix_ss( int j1, int j2, int im1, int im2,
    complex_array& in_cm, int64_t nrow, int itrp )
{
  int i1, i2, icomp, ii1, i, il, ii2, jj1, jl, jj2;
  nec_float t1xi, t1yi, t1zi, t2xi, t2yi, t2zi, xi, yi, zi;
  nec_complex g11, g12, g21, g22;
  
  i1=( im1+1)/2;
  i2=( im2+1)/2;
  icomp= i1*2-3;
  ii1 = -2;
  if ( icomp+2 < im1)
    ii1 = -3;
  
  /* loop over observation patches */
  il = -1;
  for( i = i1; i <= i2; i++ )
  {
    il++;
    icomp += 2;
    ii1 += 2;
    ii2 = ii1+1;
  
    t1xi= m_geometry->t1x[il]* m_geometry->psalp[il];
    t1yi= m_geometry->t1y[il]* m_geometry->psalp[il];
    t1zi= m_geometry->t1z[il]* m_geometry->psalp[il];
    t2xi= m_geometry->t2x[il]* m_geometry->psalp[il];
    t2yi= m_geometry->t2y[il]* m_geometry->psalp[il];
    t2zi= m_geometry->t2z[il]* m_geometry->psalp[il];
    xi= m_geometry->px[il];
    yi= m_geometry->py[il];
    zi= m_geometry->pz[il];
  
    /* loop over source patches */
    jj1 = -2;
    for(int j = j1; j <= j2; j++ )
    {
      jl=j-1;
      jj1 += 2;
      jj2 = jj1+1;
    
      m_s= m_geometry->pbi[jl];
      xj= m_geometry->px[jl];
      yj= m_geometry->py[jl];
      zj= m_geometry->pz[jl];
      t1xj= m_geometry->t1x[jl];
      t1yj= m_geometry->t1y[jl];
      t1zj= m_geometry->t1z[jl];
      t2xj= m_geometry->t2x[jl];
      t2yj= m_geometry->t2y[jl];
      t2zj= m_geometry->t2z[jl];
    
      hintg( xi, yi, zi);
    
      g11 = -( t2xi* exk+ t2yi* eyk+ t2zi* ezk);
      g12 = -( t2xi* exs+ t2yi* eys+ t2zi* ezs);
      g21 = -( t1xi* exk+ t1yi* eyk+ t1zi* ezk);
      g22 = -( t1xi* exs+ t1yi* eys+ t1zi* ezs);
    
      if ( i == j )
      {
        g11 -= .5;
        g22 += .5;
      }
    
      /* normal fill */
      if ( itrp == 0)
      {
        if ( icomp >= im1 )
        {
          in_cm[ii1+jj1*nrow]= g11;
          in_cm[ii1+jj2*nrow]= g12;
        }
      
        if ( icomp >= im2 )
          continue;
      
        in_cm[ii2+jj1*nrow]= g21;
        in_cm[ii2+jj2*nrow]= g22;
        continue;  
      } /* if ( itrp == 0) */
    
      /* transposed fill */
      if ( icomp >= im1 )
      {
        in_cm[jj1+ii1*nrow]= g11;
        in_cm[jj2+ii1*nrow]= g12;
      }
    
      if ( icomp >= im2 )
        continue;
    
      in_cm[jj1+ii2*nrow]= g21;
      in_cm[jj2+ii2*nrow]= g22;
    
    } /* for( j = j1; j <= j2; j++ ) */
  } /* for( i = i1; i <= i2; i++ ) */
}
  

/*-----------------------------------------------------------------------*/

/* computes matrix elements for e along wires due to patch current */
void nec_context::cmsw( int j1, int j2, int i1, int i2, complex_array& in_cm,
    complex_array& cw, int64_t ncw, int64_t nrow, int itrp )
{
  nec_float fsign=1.0;
  complex_array emel(9);
  
  /*! -1 offset to "m_geometry->jsno" for array indexing */
  int jsnox = m_geometry->jsno-1;
  
  if (itrp < 0)
    return;

  int k = -1;
  int icgo=0;
  
  /* observation loop */
  for (int i = i1-1; i < i2; i++ ) {
    k++;
    nec_float xi = m_geometry->x[i];
    nec_float yi= m_geometry->y[i];
    nec_float zi= m_geometry->z[i];
    nec_float cabi= m_geometry->cab[i];
    nec_float sabi= m_geometry->sab[i];
    nec_float salpi= m_geometry->salp[i];
    
    int ipch=0;
    
    if ( m_geometry->icon1[i] >= PCHCON) {
      ipch= m_geometry->icon1[i]-PCHCON;
      fsign = -1.;
    }
    
    if ( m_geometry->icon2[i] >= PCHCON) {
      ipch= m_geometry->icon2[i]-PCHCON;
      fsign=1.;
    }
    
    /* source loop */
    int jl = -1;
    for (int j = j1; j <= j2; j++ ) {
      jl += 2;
      int js = j-1;
      t1xj= m_geometry->t1x[js];
      t1yj= m_geometry->t1y[js];
      t1zj= m_geometry->t1z[js];
      t2xj= m_geometry->t2x[js];
      t2yj= m_geometry->t2y[js];
      t2zj= m_geometry->t2z[js];
      xj= m_geometry->px[js];
      yj= m_geometry->py[js];
      zj= m_geometry->pz[js];
      m_s = m_geometry->pbi[js];
    
      /* ground loop */
     // TODO weird problem to test here. replaced loop variable called ip
      for (int ipgnd = 1; ipgnd <= ground.ksymp; ipgnd++ ) {  
        if ( ((ipch == j) || (icgo != 0)) && (ipgnd != 2) ) {
          if ( icgo <= 0 ) {
            ASSERT(ipgnd != 2);
            pcint( xi, yi, zi, cabi, sabi, salpi, emel);
          
            nec_float angle = pi()* m_geometry->segment_length[i]* fsign;
            nec_float pxl = sin(angle);
            nec_float pyl = cos(angle);

            exc= emel[8]* fsign;
          
            m_geometry->trio(i+1);
          
            int il= i-ncw;
            if ( i < m_geometry->np)
              il += (il/m_geometry->np)*2*m_geometry->mp;
          
            if ( itrp == 0 )
              cw[k+il*nrow] += exc*( m_geometry->ax[jsnox]+ m_geometry->bx[jsnox]* pxl+ m_geometry->cx[jsnox]* pyl);
            else
              cw[il+k*nrow] += exc*( m_geometry->ax[jsnox]+ m_geometry->bx[jsnox]* pxl+ m_geometry->cx[jsnox]* pyl);
          } /* if ( icgo <= 0 ) */
        
          if ( itrp == 0) {
            in_cm[k+(jl-1)*nrow]= emel[icgo];
            in_cm[k+jl*nrow]    = emel[icgo+4];
          } else {
            in_cm[(jl-1)+k*nrow]= emel[icgo];
            in_cm[jl+k*nrow]    = emel[icgo+4];
          }
        
          icgo++;
          if ( icgo == 4)
            icgo=0;
        
          continue;
        } /* if ( ((ipch == (j+1)) || (icgo != 0)) && (ipgnd != 2) ) */
      
        unere( xi, yi, zi, ipgnd == 2);
      
        /* normal fill */
        if ( itrp == 0) {
          in_cm[k+(jl-1)*nrow] += exk* cabi+ eyk* sabi+ ezk* salpi;
          in_cm[k+jl*nrow]     += exs* cabi+ eys* sabi+ ezs* salpi;
        } else {
          /* transposed fill */
          in_cm[(jl-1)+k*nrow] += exk* cabi+ eyk* sabi+ ezk* salpi;
          in_cm[jl+k*nrow]     += exs* cabi+ eys* sabi+ ezs* salpi;
        }
      } /* for( ipgnd = 1; ipgnd <= ground.ksymp; ipgnd++ ) */
    } /* for( j = j1; j <= j2; j++ ) */
  } /* for( i = i1-1; i < i2; i++ ) */
}

/*-----------------------------------------------------------------------*/

/* cmws computes matrix elements for wire-surface interactions */
void nec_context::cmws( int j, int i1, int i2, complex_array& in_cm,
    int64_t nr, complex_array& cw, int64_t nw, int itrp )
{
  UNUSED(nw);
  int ipr, js=0;
  nec_float xi, yi, zi, tx, ty, tz;
  nec_complex etk, ets, etc;
  
  j--;
  m_s= m_geometry->segment_length[j];
  m_b= m_geometry->segment_radius[j];
  xj= m_geometry->x[j];
  yj= m_geometry->y[j];
  zj= m_geometry->z[j];
  cabj= m_geometry->cab[j];
  sabj= m_geometry->sab[j];
  salpj= m_geometry->salp[j];
  
  /* observation loop */
  ipr= -1;
  for(int i = i1; i <= i2; i++ ) {
    ipr++;
    int ipatch=(i+1)/2;
    int ik= i-( i/2)*2;

    if ( (ik != 0) || (ipr == 0) ) {
      js= ipatch-1;
      xi= m_geometry->px[js];
      yi= m_geometry->py[js];
      zi= m_geometry->pz[js];
      hsfld( xi, yi, zi,0.);

      if ( ik != 0 ) {
        tx= m_geometry->t2x[js];
        ty= m_geometry->t2y[js];
        tz= m_geometry->t2z[js];
      } else {
        tx= m_geometry->t1x[js];
        ty= m_geometry->t1y[js];
        tz= m_geometry->t1z[js];
      }
    } /* if ( (ik != 0) || (ipr == 0) ) */
    else {
      tx= m_geometry->t1x[js];
      ty= m_geometry->t1y[js];
      tz= m_geometry->t1z[js];
    } /* if ( (ik != 0) || (ipr == 0) ) */

    etk = -( exk* tx+ eyk* ty+ ezk* tz)* m_geometry->psalp[js];
    ets = -( exs* tx+ eys* ty+ ezs* tz)* m_geometry->psalp[js];
    etc = -( exc* tx+ eyc* ty+ ezc* tz)* m_geometry->psalp[js];

    /* fill matrix elements.  element locations */
    /* determined by connection data. */
    
    switch (itrp) {
      case 0:   /* normal fill */
        for(int ij = 0; ij < m_geometry->jsno; ij++ ) {
          int jx= m_geometry->jco[ij]-1;
          in_cm[ipr+jx*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
          /* CM(IPR,JX)=CM(IPR,JX)+ETK*AX(IJ)+ETS*BX(IJ)+ETC*CX(IJ) */
        }
        break;
      
      case 2:  /* transposed fill - c(ws) and d(ws)prime (=cw) */
        for (int ij = 0; ij < m_geometry->jsno; ij++ ) {
          int jx= m_geometry->jco[ij]-1;
          if ( jx < nr) {
            in_cm[jx+ipr*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
          } else {
            jx -= nr;
            cw[jx+ipr*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
          }
        }
        break;
        
      default:  /* transposed fill */
        for(int ij = 0; ij < m_geometry->jsno; ij++ ) {
          int jx= m_geometry->jco[ij]-1;
          in_cm[jx+ipr*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
          /* CM(JX,IPR)=CM(JX,IPR)+ETK*AX(IJ)+ETS*BX(IJ)+ETC*CX(IJ) */
        }
        break;
    }
  } /* for( i = i1; i <= i2; i++ ) */
}

/*-----------------------------------------------------------------------*/

/* cmww computes matrix elements for wire-wire interactions */
void nec_context::cmww( int j, int i1, int i2, complex_array& in_cm,
    int64_t nr, complex_array& cw, int64_t nw, int itrp) 
{
  UNUSED(nw);
  int i, jx;
  nec_float xi, yi, zi, ai, cabi, sabi, salpi;
  nec_complex etk, ets, etc;
  
  /* set source segment parameters */
  jx = j;
  j--;
  m_s= m_geometry->segment_length[j];
  m_b= m_geometry->segment_radius[j];
  xj= m_geometry->x[j];
  yj= m_geometry->y[j];
  zj= m_geometry->z[j];
  cabj= m_geometry->cab[j];
  sabj= m_geometry->sab[j];
  salpj= m_geometry->salp[j];
  
  /* Decide whether ext. t.w. approx. can be used */
  if ( m_use_exk == true) {
    int ipr = m_geometry->icon1[j];
#if 0
    if (ipr > PCHCON) ind1 = 0; // trap bug also in original fortran
    else if (ipr == 0) ind1 = 1;
    else if ((ipr == jx) && ( cabj*cabj + sabj*sabj > 1.e-8))
            ind1=2;
    else if (ipr < 0) {
      ind1 = 2;
      int iprx= -ipr - 1;
      if ( m_geometry->icon1[iprx] == -jx )
        ind1 = m_geometry->test_ek_approximation(j,iprx);
    } else {  /* if ( ipr > 0 ) */
      if ( ipr != jx ) {
        ind1=2;
        int iprx = ipr-1;
        if ( m_geometry->icon2[iprx] == jx )
          ind1 = m_geometry->test_ek_approximation(j,iprx);
      }
      else if ( (cabj* cabj+ sabj* sabj) > 1.e-8)
        ind1=2;
      else
        ind1=0;
    } /* if ( ipr < 0 ) */
#endif
    {
      int iprx = std::abs(ipr) - 1;
      ind1 = 2;
      if (ipr > PCHCON)
        ind1 = 0;
      else if (ipr == 0)
        ind1 = 1;
      else if ((ipr == jx) && (cabj*cabj + sabj*sabj <= 1.e-8))
        ind1 = 0;
      else if ((ipr > 0) && (m_geometry->icon2[iprx] == jx))
        ind1 = m_geometry->test_ek_approximation(j, ipr - 1);
      else if ((ipr < 0) && (m_geometry->icon1[iprx] == -jx))
        ind1 = m_geometry->test_ek_approximation(j, -ipr - 1);
    }

    ipr = m_geometry->icon2[j];

    if (ipr > PCHCON) ind2 = 2; // trap bug also in original fortran
    else if (ipr == 0) ind2 = 1;
    else if (ipr < 0) {
      int iprx = -ipr - 1;
      if ( -m_geometry->icon2[iprx] != jx ) {
        ind2=2;
      } else {
        ind2 = m_geometry->test_ek_approximation(j,iprx);
      }
    } else { /* if ( ipr>= 0 ) */
      if ( ipr != jx ) {
        int iprx = ipr-1;
        if ( m_geometry->icon1[iprx] != jx ) {
          ind2=2;
        } else {
          ind2 = m_geometry->test_ek_approximation(j,iprx);
        }
      }
      else if ( (cabj* cabj+ sabj* sabj) > 1.e-8) {
        ind2=2;
      } else {
        ind2=0;
      }
    } /* if ( ipr < 0 ) */
  } /* if ( m_use_exk == true) */
  
  /* observation loop */
  int ipr = -1;
  for( i = i1-1; i < i2; i++ ) {
    ipr++;
    xi= m_geometry->x[i];
    yi= m_geometry->y[i];
    zi= m_geometry->z[i];
    ai= m_geometry->segment_radius[i];
    cabi= m_geometry->cab[i];
    sabi= m_geometry->sab[i];
    salpi= m_geometry->salp[i];

    efld( xi, yi, zi, ai, i != j);

    etk= exk* cabi+ eyk* sabi+ ezk* salpi;
    ets= exs* cabi+ eys* sabi+ ezs* salpi;
    etc= exc* cabi+ eyc* sabi+ ezc* salpi;

    /* fill matrix elements. element locations */
    /* determined by connection data. */

    /* normal fill */
    if ( itrp == 0) {
      for (int ij = 0; ij < m_geometry->jsno; ij++ ) {
        jx = m_geometry->jco[ij]-1;
        in_cm[ipr+jx*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
      }
      continue;
    }

    /* transposed fill */
    if ( itrp != 2) {
      for (int ij = 0; ij < m_geometry->jsno; ij++ ) {
        jx= m_geometry->jco[ij]-1;
        in_cm[jx+ipr*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
      }
      continue;
    }

    /* trans. fill for c(ww) - test for elements for d(ww)prime.  (=cw) */
    for (int ij = 0; ij < m_geometry->jsno; ij++ ) {
      jx= m_geometry->jco[ij]-1;
      if ( jx < nr) {
        in_cm[jx+ipr*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
      } else {
        jx -= nr;
        cw[jx+ipr*nr] += etk* m_geometry->ax[ij]+ ets* m_geometry->bx[ij]+ etc* m_geometry->cx[ij];
      }
    } /* for( ij = 0; ij < m_geometry->jsno; ij++ ) */

  } /* for( i = i1-1; i < i2; i++ ) */
}


/*-----------------------------------------------------------------------*/

/* couple computes the maximum coupling between pairs of segments. */
void nec_context::couple( complex_array& in_currents, nec_float in_wavelength )
{
  int j1, j2, l1, itt1, itt2, its1, its2, isg1, isg2, npm1;
  nec_float dbc, c, gmax;
  nec_complex y11, y12, y22, yl, yin, zl, zin, rho;
  
  if ( (voltage_source_count != 1) || (nvqd != 0) )
    return;
  
  int seg_num = m_geometry->get_segment_number( nctag[icoup], ncseg[icoup]);
  if ( seg_num != source_segment_array[0] )
    return;
  
  zin= source_voltage_array[0];
  icoup++;
  y11a.resize(icoup);
  y11a[icoup-1]= in_currents[seg_num-1]*in_wavelength/zin;
  
  l1=(icoup-1)*(ncoup-1);
  for (int i = 0; i < ncoup; i++ )
  {
    if ( (i+1) == icoup)
      continue;
  
    l1++;
    y12a.resize(l1);
    int k = m_geometry->get_segment_number( nctag[i], ncseg[i]);
    y12a[l1-1]= in_currents[k-1]* in_wavelength/ zin;
  }
  
  if ( icoup < ncoup)
    return;
  
  m_output.endl(3);
  m_output.nec_printf(
    "                        -----------"
    " ISOLATION DATA -----------\n\n"
    " ------- COUPLING BETWEEN ------     MAXIMUM    "
    " ---------- FOR MAXIMUM COUPLING ----------\n"
    "            SEG              SEG    COUPLING  LOAD"
    " IMPEDANCE (2ND SEG)         INPUT IMPEDANCE \n"
    " TAG  SEG   No:   TAG  SEG   No:      (DB)       "
    " REAL     IMAGINARY         REAL       IMAGINARY" );
  
  npm1= ncoup-1;
  
  for (int i = 0; i < npm1; i++ )
  {
    itt1= nctag[i];
    its1= ncseg[i];
    isg1= m_geometry->get_segment_number( itt1, its1);
    l1= i+1;
  
    for (int j = l1; j < ncoup; j++ )
    {
      itt2= nctag[j];
      its2= ncseg[j];
      isg2= m_geometry->get_segment_number( itt2, its2);
      j1= j+ i* npm1-1;
      j2= i+ j* npm1;
      y11= y11a[i];
      y22= y11a[j];
      y12=.5*( y12a[j1]+ y12a[j2]);
      yin= y12* y12;
      dbc= abs( yin);
      c= dbc/(2.* real( y11)* real( y22)- real( yin));
    
      if ( (c >= 0.0) && (c <= 1.0) )
      {
        if ( c >= .01 )
          gmax=(1.- sqrt(1.- c*c))/c;
        else
          gmax=.5*( c+.25* c* c* c);
      
        rho= gmax* conj( yin)/ dbc;
        yl=((1.- rho)/(1.+ rho)+1.)* real( y22)- y22;
        zl=1./ yl;
        yin= y11- yin/( y22+ yl);
        zin=1./ yin;
        dbc= db10( gmax);
      
        m_output.endl();
        m_output.nec_printf(
          " %4d %4d %5d  %4d %4d %5d  %9.3f"
          "  %12.5E %12.5E  %12.5E %12.5E",
          itt1, its1, isg1, itt2, its2, isg2, dbc,
          real(zl), imag(zl), real(zin), imag(zin) );
      
        continue;
      
      } /* if ( (c >= 0.0) && (c <= 1.0) ) */
    
      m_output.endl();
      m_output.nec_printf(
        " %4d %4d %5d   %4d %4d %5d  **ERROR** "
        "COUPLING IS NOT BETWEEN 0 AND 1. (= %12.5E)",
        itt1, its1, isg1, itt2, its2, isg2, c );
    
    } /* for( j = l1; j < ncoup; j++ ) */
  
  } /* for( i = 0; i < npm1; i++ ) */
}


/*! \brief Compute near E-fields of a segment due to constant sine and cosine currents distributions.  The ground effect is included.

\param xi,yi,zi x,y,z components of the field evaluation point.
\param ai  Radius of the segment on which the field is evaluated.
\param on_source_segment  Flag to indicate if the field evaluation point is not on the source segment (i != j)
*/
void nec_context::efld( nec_float xi, nec_float yi, nec_float zi, nec_float ai, bool not_on_source_segment)
{
  #define  txk  egnd[0]
  #define  tyk  egnd[1]
  #define  tzk  egnd[2]
  #define  txs  egnd[3]
  #define  tys  egnd[4]
  #define  tzs  egnd[5]
  #define  txc  egnd[6]
  #define  tyc  egnd[7]
  #define  tzc  egnd[8]
  
  nec_float salpr, zij, zp, rhox;
  nec_float rhoy, rhoz, rh, r, rmag, cth, px, py;
  nec_float xymag, xspec, yspec, rhospc, dmin, shaf;
  nec_complex epx, epy, refs, refps, zrsin, zratx, zscrn;
  nec_complex tezs, ters, tezc, terc, tezk, terk;
  complex_array egnd(9);
  
  nec_float xij = xi - xj;
  nec_float yij = yi - yj;
  bool ijx = not_on_source_segment;
  
  {
    // calculate the direct field.
    salpr= salpj;
    zij= zi- zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;
  
    rh = sqrt( rhox*rhox + rhoy*rhoy + rhoz*rhoz + ai*ai);
    if ( rh <= 1.e-10)
    {
      rhox=0.0;
      rhoy=0.0;
      rhoz=0.0;
    }
    else
    {
      rhox= rhox/ rh;
      rhoy= rhoy/ rh;
      rhoz= rhoz/ rh;
    }
  
    /* lumped current element approx. for large separations */
    r= sqrt( zp*zp + rh*rh);
    if ( r >= rkh)
    {
      rmag = two_pi() * r;
      cth = zp/ r;
      px = rh/ r;
      txk = nec_complex( cos( rmag),- sin( rmag));
      py= two_pi() * r* r;
      tyk= em::impedance() * cth* txk* nec_complex(1.0,-1.0/ rmag)/ py;
      tzk= em::impedance() * px* txk* nec_complex(1.0, rmag-1.0/ rmag)/(2.* py);
      tezk= tyk* cth- tzk* px;
      terk= tyk* px+ tzk* cth;
      rmag= sin( pi()* m_s)/ pi();
      tezc= tezk* rmag;
      terc= terk* rmag;
      tezk= tezk* m_s;
      terk= terk* m_s;
      txs=cplx_00();
      tys=cplx_00();
      tzs=cplx_00();
    
    } /* if ( r >= rkh) */
  
    if ( r < rkh)
    {
      /* eksc for thin wire approx. or ekscx for extended t.w. approx. */
      if ( m_use_exk == false)
        eksc( m_s, zp, rh, two_pi(), ijx, &tezs, &ters, &tezc, &terc, &tezk, &terk );
      else
        ekscx( m_b, m_s, zp, rh, two_pi(), ijx, ind1, ind2, &tezs, &ters, &tezc, &terc, &tezk, &terk);
    
      txs= tezs* cabj+ ters* rhox;
      tys= tezs* sabj+ ters* rhoy;
      tzs= tezs* salpr+ ters* rhoz;
    }
  
    txk= tezk* cabj+ terk* rhox;
    tyk= tezk* sabj+ terk* rhoy;
    tzk= tezk* salpr+ terk* rhoz;
    txc= tezc* cabj+ terc* rhox;
    tyc= tezc* sabj+ terc* rhoy;
    tzc= tezc* salpr+ terc* rhoz;
  
    exk= txk;
    eyk= tyk;
    ezk= tzk;
    exs= txs;
    eys= tys;
    ezs= tzs;
    exc= txc;
    eyc= tyc;
    ezc= tzc;
  }
  
  if (ground.present())
  {
    // Now do the reflected field...
    ijx=1;
    salpr= -salpj;
    zij= zi + zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;
  
    rh = sqrt( rhox*rhox + rhoy*rhoy + rhoz*rhoz + ai*ai);
    if ( rh <= 1.e-10)
    {
      rhox=0.0;
      rhoy=0.0;
      rhoz=0.0;
    }
    else
    {
      rhox= rhox/ rh;
      rhoy= rhoy/ rh;
      rhoz= rhoz/ rh;
    }
  
    /* lumped current element approx. for large separations */
    r= sqrt( zp*zp + rh*rh);
    if ( r >= rkh)
    {
      rmag = two_pi() * r;
      cth = zp/ r;
      px = rh/ r;
      txk = nec_complex( cos( rmag),- sin( rmag));
      py= two_pi() * r* r;
      tyk= em::impedance() * cth* txk* nec_complex(1.0,-1.0/ rmag)/ py;
      tzk= em::impedance() * px* txk* nec_complex(1.0, rmag-1.0/ rmag)/(2.* py);
      tezk= tyk* cth- tzk* px;
      terk= tyk* px+ tzk* cth;
      rmag= sin( pi()* m_s)/ pi();
      tezc= tezk* rmag;
      terc= terk* rmag;
      tezk= tezk* m_s;
      terk= terk* m_s;
      txs=cplx_00();
      tys=cplx_00();
      tzs=cplx_00();
    
    } /* if ( r >= rkh) */
  
    if ( r < rkh)
    {
      /* eksc for thin wire approx. or ekscx for extended t.w. approx. */
      if ( m_use_exk == false)
        eksc( m_s, zp, rh, two_pi(), ijx, &tezs, &ters, &tezc, &terc, &tezk, &terk );
      else
        ekscx( m_b, m_s, zp, rh, two_pi(), ijx, ind1, ind2, &tezs, &ters, &tezc, &terc, &tezk, &terk);
    
      txs= tezs* cabj+ ters* rhox;
      tys= tezs* sabj+ ters* rhoy;
      tzs= tezs* salpr+ ters* rhoz;
    }
  
    txk= tezk* cabj+ terk* rhox;
    tyk= tezk* sabj+ terk* rhoy;
    tzk= tezk* salpr+ terk* rhoz;
    txc= tezc* cabj+ terc* rhox;
    tyc= tezc* sabj+ terc* rhoy;
    tzc= tezc* salpr+ terc* rhoz;
  
    ASSERT(ground.is_valid());
    if (ground.type_finite_reflection())
    {
      zratx= ground.zrati;
      rmag= r;
      xymag= sqrt( xij* xij+ yij* yij);
    
      /* set parameters for radial wire ground screen. */
      if (  ground.get_radial_wire_count() != 0)
      {
        xspec=( xi* zj+ zi* xj)/( zi+ zj);
        yspec=( yi* zj+ zi* yj)/( zi+ zj);
        rhospc= sqrt( xspec*xspec + yspec*yspec + ground.t2*ground.t2);
      
        if ( rhospc <= ground.get_radial_wire_length_wavelengths())
        {
          zscrn= ground.m_t1* rhospc* log( rhospc/ ground.t2);
          zratx=( zscrn* ground.zrati)/( em::impedance() * ground.zrati+ zscrn);
        }
      } /* if (  ground.get_radial_wire_count() != 0) */
    
      /* Calculation of reflection coefficients when ground is specified. */
      if ( xymag <= 1.0e-6)
      {
        px=0.;
        py=0.;
        cth=1.;
        zrsin=cplx_10();
      }
      else
      {
        px = - yij/ xymag;
        py= xij/ xymag;
        cth= zij/ rmag;
        zrsin= sqrt(1.0 - zratx*zratx*(1.0 - cth*cth) );
      } /* if ( xymag <= 1.0e-6) */
    
      refs=( cth- zratx* zrsin)/( cth+ zratx* zrsin);
      refps = -( zratx* cth- zrsin)/( zratx* cth+ zrsin);
      refps= refps- refs;
      epy= px* txk+ py* tyk;
      epx= px* epy;
      epy= py* epy;
      txk= refs* txk+ refps* epx;
      tyk= refs* tyk+ refps* epy;
      tzk= refs* tzk;
      epy= px* txs+ py* tys;
      epx= px* epy;
      epy= py* epy;
      txs= refs* txs+ refps* epx;
      tys= refs* tys+ refps* epy;
      tzs= refs* tzs;
      epy= px* txc+ py* tyc;
      epx= px* epy;
      epy= py* epy;
      txc= refs* txc+ refps* epx;
      tyc= refs* tyc+ refps* epy;
      tzc= refs* tzc;
    
    } /* if (ground.type_finite_reflection()) */
    
    exk -= txk* ground.frati;
    eyk -= tyk* ground.frati;
    ezk -= tzk* ground.frati;
    exs -= txs* ground.frati;
    eys -= tys* ground.frati;
    ezs -= tzs* ground.frati;
    exc -= txc* ground.frati;
    eyc -= tyc* ground.frati;
    ezc -= tzc* ground.frati;
  }
  
  if (false == ground.type_sommerfeld_norton())
    return;
  
  /* field due to ground using Sommerfeld/Norton */
  sn = norm(cabj, sabj);
  if ( sn >= 1.0e-5) {
    xsn= cabj/ sn;
    ysn= sabj/ sn;
  } else {
    sn=0.0;
    xsn=1.0;
    ysn=0.0;
  }
  
  /* displace observation point for thin wire approximation */
  zij= zi+ zj;
  salpr = - salpj;
  rhox= sabj* zij- salpr* yij;
  rhoy= salpr* xij- cabj* zij;
  rhoz= cabj* yij- sabj* xij;
  rh = norm(rhox, rhoy, rhoz); // rhox* rhox+ rhoy* rhoy+ rhoz* rhoz;
  
  if ( rh <= 1.e-10) {
    xo= xi- ai* ysn;
    yo= yi+ ai* xsn;
    zo= zi;
  } else {
    rh= ai/ sqrt( rh);
    if ( rhoz < 0.0)
      rh = - rh;
    xo= xi+ rh* rhox;
    yo= yi+ rh* rhoy;
    zo= zi+ rh* rhoz;
  } /* if ( rh <= 1.e-10) */
  
  r = xij*xij + yij*yij + zij*zij;
  
  if ( r <= .95) {
    /* Field from interpolation is integrated over segment */
    isnor=1;
    dmin = norm(exk) + norm(eyk) + norm(ezk);
  
    dmin = 0.01* sqrt(dmin);
    shaf = 0.5* m_s;
    rom2(-shaf, shaf, egnd, dmin);
  } else {
    /* Norton field equations and lumped current element approximation */
    isnor=2;
    sflds(0., egnd);
  } /* if ( r <= .95) */
  
  if ( r > .95) {
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rh= r- zp* zp;
    if ( rh <= 1.e-10)
      dmin = 0.0;
    else
      dmin = sqrt( rh/( rh+ ai* ai));
  
    if ( dmin <= .95) {
      px=1.- dmin;
      terk=( txk* cabj+ tyk* sabj+ tzk* salpr)* px;
      txk= dmin* txk+ terk* cabj;
      tyk= dmin* tyk+ terk* sabj;
      tzk= dmin* tzk+ terk* salpr;
      ters=( txs* cabj+ tys* sabj+ tzs* salpr)* px;
      txs= dmin* txs+ ters* cabj;
      tys= dmin* tys+ ters* sabj;
      tzs= dmin* tzs+ ters* salpr;
      terc=( txc* cabj+ tyc* sabj+ tzc* salpr)* px;
      txc= dmin* txc+ terc* cabj;
      tyc= dmin* tyc+ terc* sabj;
      tzc= dmin* tzc+ terc* salpr;
    }
  } /* if ( r > .95) */
  
  exk= exk+ txk;
  eyk= eyk+ tyk;
  ezk= ezk+ tzk;
  exs= exs+ txs;
  eys= eys+ tys;
  ezs= ezs+ tzs;
  exc= exc+ txc;
  eyc= eyc+ tyc;
  ezc= ezc+ tzc;
}

/*-----------------------------------------------------------------------*/

/*
  Compute E-field of sine, cosine, and constant
  current filaments by thin wire approximation.
*/
void nec_context::eksc( nec_float s, nec_float z, nec_float rh, nec_float xk, int ij,
    nec_complex *in_ezs, nec_complex *ers, nec_complex *in_ezc,
    nec_complex *erc, nec_complex *in_ezk, nec_complex *erk )
{
  static nec_complex __const1(0.0,4.771341189);
  
  // set some global variables for the gf() method. These should be moved to parameters of gf (via intx)
  ija = ij;
  zpk = xk* z;
  
  nec_float rhk = xk * rh;
  
  rkb2 = rhk * rhk;
  
  nec_float sh =.5* s;
  nec_float shk = xk * sh;
  nec_float ss = sin(shk);
  nec_float cs = cos(shk);
  nec_float z2a = sh - z;
  nec_float z1a = -(sh + z);
  
  nec_complex gz1, gz2, gp1, gp2;
  gx( z1a, rh, xk, &gz1, &gp1);
  gx( z2a, rh, xk, &gz2, &gp2);  
  nec_complex gzp1 = gp1 * z1a;
  nec_complex gzp2 = gp2 * z2a;
  
  *in_ezs=  __const1*(( gz2- gz1)* cs* xk-( gzp2+ gzp1)* ss);
  *in_ezc = - __const1*(( gz2+ gz1)* ss* xk+( gzp2- gzp1)* cs);
  *erk= __const1*( gp2- gp1)* rh;
  
  nec_float cint, sint;
  intx(-shk, shk, rhk, ij, &cint, &sint);
  *in_ezk = - __const1*( gzp2- gzp1+ xk* xk* nec_complex( cint,- sint));
  
  if ( rh >= 1.0e-10)
  {
    gzp1 = gzp1 * z1a;
    gzp2 = gzp2 * z2a;
    
    *ers = - __const1*(( gzp2+ gzp1+ gz2+ gz1)*
      ss-( z2a* gz2- z1a* gz1)* cs*xk)/ rh;
    *erc = - __const1*(( gzp2- gzp1+ gz2- gz1)*
      cs+( z2a* gz2+ z1a* gz1)* ss*xk)/ rh;
    return;
  }
  
  *ers = cplx_00();
  *erc = cplx_00();
}

/*-----------------------------------------------------------------------*/

/*
  Compute e field of sine, cosine, and constant current
  filaments by extended thin wire approximation.
*/
void nec_context::ekscx( nec_float bx, nec_float s, nec_float z,
    nec_float rhx, nec_float xk, int ij, int inx1, int inx2,
    nec_complex *in_ezs, nec_complex *ers, nec_complex *in_ezc,
    nec_complex *erc, nec_complex *in_ezk, nec_complex *erk )
{
  static nec_complex __const1(0.0,4.771341189);
  int ira;
  nec_float b, rh, sh, rhk, shk, ss, cs, z1a;
  nec_float z2a, a2, bk, bk2, cint, sint;
  nec_complex gz1, gz2, gzp1, gzp2, gr1, gr2;
  nec_complex grp1, grp2, grk1, grk2, gzz1, gzz2;
  
  if ( rhx >= bx)
  {
    rh= rhx;
    b= bx;
    ira=0;
  }
  else
  {
    rh= bx;
    b= rhx;
    ira=1;
  }
  
  sh=.5* s;
  
  // set some global variables for the gf() method. These should be moved to parameters of gf (via intx)
  ija= ij;
  zpk= xk* z;
  rhk= xk* rh;
  rkb2= rhk* rhk;
  
  shk= xk* sh;
  ss= sin( shk);
  cs= cos( shk);
  z2a= sh- z;
  z1a = -( sh+ z);
  a2= b* b;
  
  if ( inx1 != 2)
    gxx( z1a, rh, b, a2, xk, ira, &gz1, &gzp1, &gr1, &grp1, &grk1, &gzz1);
  else
  {
    gx( z1a, rhx, xk, &gz1, &grk1);
    gzp1= grk1* z1a;
    gr1= gz1/ rhx;
    grp1= gzp1/ rhx;
    grk1= grk1* rhx;
    gzz1= cplx_00();
  }
  
  if ( inx2 != 2)
    gxx( z2a, rh, b, a2, xk, ira, &gz2, &gzp2, &gr2, &grp2, &grk2, &gzz2);
  else
  {
    gx( z2a, rhx, xk, &gz2, &grk2);
    gzp2= grk2* z2a;
    gr2= gz2/ rhx;
    grp2= gzp2/ rhx;
    grk2= grk2* rhx;
    gzz2= cplx_00();
  }
  
  *in_ezs= __const1*(( gz2- gz1)* cs* xk-( gzp2+ gzp1)* ss);
  *in_ezc = - __const1*(( gz2+ gz1)* ss* xk+( gzp2- gzp1)* cs);
  *ers = - __const1*(( z2a* grp2+ z1a* grp1+ gr2+ gr1)*ss
    -( z2a* gr2- z1a* gr1)* cs* xk);
  *erc = - __const1*(( z2a* grp2- z1a* grp1+ gr2- gr1)*cs
    +( z2a* gr2+ z1a* gr1)* ss* xk);
  *erk= __const1*( grk2- grk1);
  intx(- shk, shk, rhk, ij, &cint, &sint);
  bk= b* xk;
  bk2= bk* bk*.25;
  *in_ezk = - __const1*( gzp2- gzp1+ xk* xk*(1.- bk2)*
    nec_complex( cint,- sint)-bk2*( gzz2- gzz1));
}

/*-----------------------------------------------------------------------*/

/*!\brief Fills the array e with the negative of the electric field tangent to the segments and the tangential magnetic field on the surfaces

\param e is the right hand side of the matrix equation.
*/
void nec_context::etmns( nec_float p1, nec_float p2, nec_float p3, nec_float p4,
    nec_float p5, nec_float p6, nec_float incident_amplitude, enum excitation_type excite_type, complex_array& e )
{
  nec_float cth, sth, cph, sph, cet, set, pxl, pyl, pzl, wx;
  nec_float wy, wz, qx, qy, qz, ds, dsh, rs, r;
  nec_complex er, et, ezh, erh, rrv, rrh, tt1, tt2;
  
  int n = m_geometry->n_segments;
  int m = m_geometry->m;
  
//  int neq= n+2*m;
  ASSERT(neq == n+2*m);
  
  nqds=0;
  
  ASSERT(excite_type >= 0);
    
  /* applied field of voltage sources for transmitting case */
  if ( (excite_type == EXCITATION_VOLTAGE) || (excite_type == EXCITATION_VOLTAGE_DISC) )
  {
    vector_fill(e,0,neq,cplx_00());
  
    for (int i = 0; i < voltage_source_count; i++ )
    {
      int source_index = source_segment_array[i]-1;
      e[source_index] = -source_voltage_array[i]/(m_geometry->segment_length[source_index]* _wavelength);
    }
  
    if ( nvqd == 0)
      return;
  
    for (int i = 0; i < nvqd; i++ )
    {
      qdsrc( ivqd[i], vqd[i], e);
    }
    return;
  } /* if ( (excite_type == EXCITATION_VOLTAGE) || (excite_type == EXCITATION_VOLTAGE_DISC) ) */

  /* incident plane wave, linearly polarized. */
  if ((excite_type == EXCITATION_LINEAR) ||
    (excite_type == EXCITATION_CIRC_RIGHT) ||
    (excite_type == EXCITATION_CIRC_LEFT))
  {
    cth= cos( p1);
    sth= sin( p1);
    cph= cos( p2);
    sph= sin( p2);
    cet= cos( p3);
    set= sin( p3);
    pxl= cth* cph* cet- sph* set;
    pyl= cth* sph* cet+ cph* set;
    pzl = - sth* cet;
    wx = - sth* cph;
    wy = - sth* sph;
    wz = - cth;
    qx= wy* pzl- wz* pyl;
    qy= wz* pxl- wx* pzl;
    qz= wx* pyl- wy* pxl;

    if (ground.present())
    {
      if (ground.type_perfect())
      {
        rrv = -cplx_10();
        rrh = -cplx_10();
      }
      else
      {
        rrv= sqrt(1.0 - ground.get_zrati_sqr() * sth*sth);
        rrh= ground.zrati* cth;
        rrh=( rrh- rrv)/( rrh+ rrv);
        rrv= ground.zrati* rrv;
        rrv = -( cth- rrv)/( cth+ rrv);
      }
    }
  
    if ( excite_type == EXCITATION_LINEAR)
    {
      if ( n != 0)
      {
        for (int i = 0; i < n; i++ )
        {
          nec_float arg = - two_pi() *( wx* m_geometry->x[i]+ wy* m_geometry->y[i]+ wz* m_geometry->z[i]);
          nec_complex e_amplitude(cos(arg), sin(arg));
          
          e[i] = -(pxl*m_geometry->cab[i]+ pyl*m_geometry->sab[i]+ pzl*m_geometry->salp[i]) * e_amplitude * incident_amplitude;
        }
      
        if (ground.present())
        {
          tt1=( pyl* cph- pxl* sph)*( rrh- rrv);
          nec_complex cx= rrv* pxl- tt1* sph;
          nec_complex cy= rrv* pyl+ tt1* cph;
          nec_complex cz = - rrv* pzl;
        
          for (int i = 0; i < n; i++ )
          {
            nec_float arg = - two_pi()*( wx* m_geometry->x[i]+ wy* m_geometry->y[i]- wz* m_geometry->z[i]);
            nec_complex e_amplitude(cos(arg), sin(arg));
            e[i] -= ( cx* m_geometry->cab[i]+ cy* m_geometry->sab[i]+ cz* m_geometry->salp[i])* e_amplitude * incident_amplitude;
          }
        }
      } /* if ( n != 0) */
  
      if ( m == 0) // if no surface patches we're done!
        return;
  
      {
        int i= -1;
        int i1= n-2;
        for (int is = 0; is < m; is++ )
        {
          i++;
          ASSERT(is == i);
          i1 += 2;
          int i2 = i1+1;
          nec_float arg = -m_geometry->patch_angle(i,wx,wy,wz);
          ASSERT( arg == (- two_pi()*( wx* m_geometry->px[i]+ wy* m_geometry->py[i]+ wz* m_geometry->pz[i])));
          
          nec_complex e_amplitude = nec_complex(cos(arg), sin(arg)) * m_geometry->psalp[i] * em::inverse_impedance();
          e[i2]=( qx* m_geometry->t1x[i]+ qy* m_geometry->t1y[i]+ qz* m_geometry->t1z[i])* e_amplitude * incident_amplitude;
          e[i1]=( qx* m_geometry->t2x[i]+ qy* m_geometry->t2y[i]+ qz* m_geometry->t2z[i])* e_amplitude * incident_amplitude;
        }
      }
      
      if (ground.present())
      {
        tt1=( qy* cph- qx* sph)*( rrv- rrh);
        nec_complex cx = -( rrh* qx- tt1* sph);
        nec_complex cy = -( rrh* qy+ tt1* cph);
        nec_complex cz= rrh* qz;
        
        int i= -1;
        int i1= n-2;
        for (int is = 0; is < m; is++ )
        {
          i++;
          i1 += 2;
          int i2 = i1+1;
          nec_float arg = -m_geometry->patch_angle(i,wx,wy,wz);
          ASSERT((-two_pi()*( wx* m_geometry->px[i]+ wy* m_geometry->py[i]- wz* m_geometry->pz[i])) == arg);
          
          nec_complex e_amplitude = cplx_exp(arg) * m_geometry->psalp[i] * em::inverse_impedance();
          
          e[i2]= e[i2]+( cx* m_geometry->t1x[i]+ cy* m_geometry->t1y[i]+ cz* m_geometry->t1z[i])* e_amplitude * incident_amplitude;
          e[i1]= e[i1]+( cx* m_geometry->t2x[i]+ cy* m_geometry->t2y[i]+ cz* m_geometry->t2z[i])* e_amplitude * incident_amplitude;
        }
      }
      return;
    } /* if ( excite_type == EXCITATION_LINEAR) */

    /* incident plane wave, elliptic polarization. */
    tt1 = -(cplx_01())* p6;
    if ( excite_type == EXCITATION_CIRC_LEFT)
      tt1 = - tt1;
  
    if ( n != 0)
    {
      nec_complex cx= pxl+ tt1* qx;
      nec_complex cy= pyl+ tt1* qy;
      nec_complex cz= pzl+ tt1* qz;
    
      for (int i = 0; i < n; i++ )
      {
        nec_float arg = - two_pi()*( wx* m_geometry->x[i]+ wy* m_geometry->y[i]+ wz* m_geometry->z[i]);
        nec_complex e_amplitude(cos(arg), sin(arg));
        e[i] = -(cx* m_geometry->cab[i] + cy* m_geometry->sab[i] + cz*m_geometry->salp[i])* e_amplitude * incident_amplitude;
      }
      
      
      if (ground.present())
      {
        tt2=( cy* cph- cx* sph)*( rrh- rrv);
        nec_complex ccx= rrv* cx- tt2* sph;
        nec_complex ccy= rrv* cy+ tt2* cph;
        nec_complex ccz = - rrv* cz;
      
        for (int i = 0; i < n; i++ )
        {
          nec_float arg = - two_pi()*( wx* m_geometry->x[i]+ wy* m_geometry->y[i]- wz* m_geometry->z[i]);
          nec_complex e_amplitude(cos(arg), sin(arg));
          e[i] -= (ccx* m_geometry->cab[i]+ ccy* m_geometry->sab[i]+ ccz* m_geometry->salp[i]) * e_amplitude * incident_amplitude;
        }
      }
    } /* if ( n != 0) */

    if ( m == 0)
      return;

    nec_complex cx= qx- tt1* pxl;
    nec_complex cy= qy- tt1* pyl;
    nec_complex cz= qz- tt1* pzl;
  
    {
      int i= -1;
      int i1= n-2;
      for (int is = 0; is < m; is++ )
      {
        i++;
        i1 += 2;
        int i2 = i1+1;
        
        nec_float arg = -m_geometry->patch_angle(i,wx,wy,wz);
        ASSERT(arg == -two_pi()*( wx* m_geometry->px[i]+ wy* m_geometry->py[i]+ wz* m_geometry->pz[i]));

        nec_complex e_amplitude = nec_complex(cos(arg), sin(arg)) * m_geometry->psalp[i] * em::inverse_impedance();
        
        e[i2]=( cx* m_geometry->t1x[i]+ cy* m_geometry->t1y[i]+ cz* m_geometry->t1z[i])* e_amplitude * incident_amplitude;
        e[i1]=( cx* m_geometry->t2x[i]+ cy* m_geometry->t2y[i]+ cz* m_geometry->t2z[i])* e_amplitude * incident_amplitude;
      }
    }
    
    if (ground.present())
    {
      tt1=( cy* cph- cx* sph)*( rrv- rrh);
      cx = -( rrh* cx- tt1* sph);
      cy = -( rrh* cy+ tt1* cph);
      cz= rrh* cz;

      int i= -1;
      int i1= n-2;
      for (int is=0; is < m; is++ )
      {
        i++;
        i1 += 2;
        int i2 = i1+1;
        nec_float arg = -m_geometry->patch_angle(i,wx,wy,wz);
        ASSERT(arg == -two_pi()*( wx* m_geometry->px[i]+ wy* m_geometry->py[i]- wz* m_geometry->pz[i]));
        
        nec_complex e_amplitude = nec_complex(cos(arg), sin(arg)) * m_geometry->psalp[i] * em::inverse_impedance();
        
        e[i2] += ( cx* m_geometry->t1x[i]+ cy* m_geometry->t1y[i]+ cz* m_geometry->t1z[i])* e_amplitude * incident_amplitude;
        e[i1] += ( cx* m_geometry->t2x[i]+ cy* m_geometry->t2y[i]+ cz* m_geometry->t2z[i])* e_amplitude * incident_amplitude;
      }
    }
  
    return;
  
  } /* if ( excite_type <= 3) */

  ASSERT(excite_type == EXCITATION_CURRENT);
  /* incident field of an elementary current source. */

  nec_float cosp4 = cos(p4);
  
  wx = cosp4 * cos(p5);
  wy = cosp4 * sin(p5);
  wz = sin(p4);
  ds = p6*em::impedance_over_2pi();
  dsh= p6/(2.0 * two_pi());
    
  // TODO: Split the following loop up into two loops i=1:n and i=1:m
  // to loop over the patches and the segments seperately
  
  for (int i = 0; i < m_geometry->n_plus_m(); i++ )
  {
    if ( i < n )
    {
      pxl = m_geometry->x[i] - p1;
      pyl = m_geometry->y[i] - p2;
      pzl = m_geometry->z[i] - p3;
    }
    else  // patches
    {
      int patch_index = i - n;
      pxl = m_geometry->px[patch_index] - p1;
      pyl = m_geometry->py[patch_index] - p2;
      pzl = m_geometry->pz[patch_index] - p3;
    }
    
    rs = norm2(pxl,pyl,pzl);
    if ( rs < 1.0e-30)
      continue;
  
    r = sqrt(rs);
    pxl = pxl/r;
    pyl = pyl/r;
    pzl = pzl/r;
    cth = pxl*wx + pyl*wy + pzl*wz;
    sth = sqrt(1.0 - cth*cth);
    qx = pxl - wx*cth;
    qy = pyl - wy*cth;
    qz = pzl - wz*cth;
  
    nec_float arg = norm(qx,qy,qz);
    if ( arg >= 1.e-30)
    {
      qx = qx / arg;
      qy = qy / arg;
      qz = qz / arg;
    }
    else
    {
      qx = 1.0;
      qy = 0.0;
      qz = 0.0;
    } /* if ( arg >= 1.e-30) */
  
    arg = two_pi() * r;
    tt1 = nec_complex(cos(arg), -sin(arg));
  
    if ( i < n )
    {
      tt2 = nec_complex(1.0,-1.0/arg) / rs;
      er = ds* tt1* tt2* cth;
      et = 0.5 * ds * tt1 *((cplx_01()) * two_pi()/r + tt2)* sth;
      ezh= er*cth - et*sth;
      erh= er*sth + et*cth;
      nec_complex cx = ezh*wx + erh*qx;
      nec_complex cy = ezh*wy + erh*qy;
      nec_complex cz = ezh*wz + erh*qz;
      e[i] = -(cx* m_geometry->cab[i]+ cy* m_geometry->sab[i]+ cz* m_geometry->salp[i]);
    }
    else // patches
    {
      int patch_index = i - n;
      int i1 = i + patch_index*2; // was i1 += 2; and starts at n-2
      int i2 = i1+1;
      
      pxl = wy*qz - wz*qy; // cross product here...
      pyl = wz*qx - wx*qz;
      pzl = wx*qy - wy*qx;
      
      tt2= dsh* tt1* nec_complex(1.0/r, two_pi())/ r* sth* m_geometry->psalp[patch_index];
      nec_complex cx= tt2* pxl;
      nec_complex cy= tt2* pyl;
      nec_complex cz= tt2* pzl;
      e[i2] = cx*m_geometry->t1x[patch_index] + cy*m_geometry->t1y[patch_index] + cz*m_geometry->t1z[patch_index];
      e[i1] = cx*m_geometry->t2x[patch_index] + cy*m_geometry->t2y[patch_index] + cz*m_geometry->t2z[patch_index];
    } /* if ( i < n) */
  } /* for( i = 0; i < m_geometry->n_plus_m(); i++ ) */
}


/*!\brief Computes the integrand exp(jkr)/(kr) for numerical integration.
  
  TODO: Place the ija, rkb2, zpk = xk*z  global variables as parameters (possibly also the rk parameter).
*/

void nec_context::gf( nec_float zk, nec_float *co, nec_float *si )
{
  static nec_float _gf_const0 = -1.38888889e-3;
  static nec_float _gf_const1 = 4.16666667e-2;
  
  nec_float zdk = zk - zpk;
  nec_float rk = sqrt(rkb2 + zdk*zdk);
  *si= sin(rk)/ rk;
  
  if ( ija != 0 )
  {
    *co = cos(rk)/ rk;
    return;
  }
  
  if ( rk >= .2)
  {
    *co = (cos(rk) - 1.0)/rk;
    return;
  }
  
  nec_float rks = rk*rk;
  *co = ((_gf_const0 * rks + _gf_const1)* rks - 0.5)* rk;
}


/*-----------------------------------------------------------------------*/

/* integrand for h field of a wire */
void nec_context::gh( nec_float zk, nec_float *hr, nec_float *hi)
{
  nec_float rs, r, ckr, skr, rr2, rr3;
  
  rs= zk- zpka;
  rs= rhks+ rs* rs;
  r= sqrt( rs);
  ckr= cos( r);
  skr= sin( r);
  rr2=1./ rs;
  rr3= rr2/ r;
  *hr= skr* rr2+ ckr* rr3;
  *hi= ckr* rr2- skr* rr3;
}

/*-----------------------------------------------------------------------*/

/*
  Segment end contributions for thin wire approx.
  
  This function is called a lot. We should try to 
  optimise this if at all possible.
*/
void nec_context::gx( nec_float zz, nec_float rh, nec_float xk,
    nec_complex *gz, nec_complex *gzp)
{
  nec_float r2 = zz*zz+ rh*rh;
  nec_float r = sqrt(r2);
  nec_float rkz = xk * r;
  
  nec_complex temp(cos(rkz),-sin(rkz));
  temp /= r;
  *gz = temp;
  *gzp = -nec_complex(1.0, rkz) * temp / r2;
}


/*-----------------------------------------------------------------------*/

/*! \brief Segment end contributions for extended thin wire approx.
*/
void nec_context::gxx( nec_float zz, nec_float rh, nec_float a, nec_float a2, nec_float xk, int ira,
    nec_complex *g1, nec_complex *g1p, nec_complex *g2,
    nec_complex *g2p, nec_complex *g3, nec_complex *gzp )
{
  nec_float r, r2, r4, rk, rk2, rh2, t1, t2;
  nec_complex  gz, c1, c2, c3;
  
  r2= zz* zz+ rh* rh;
  r= sqrt( r2);
  r4= r2* r2;
  rk= xk* r;
  rk2= rk* rk;
  rh2= rh* rh;
  t1=.25* a2* rh2/ r4;
  t2=.5* a2/ r2;
  c1= nec_complex(1.0, rk);
  c2=3.* c1- rk2;
  c3= nec_complex(6.0, rk)* rk2-15.* c1;
  gz= nec_complex( cos( rk),- sin( rk))/ r;
  *g2= gz*(1.+ t1* c2);
  *g1= *g2- t2* c1* gz;
  gz= gz/ r2;
  *g2p= gz*( t1* c3- c1);
  *gzp= t2* c2* gz;
  *g3= *g2p+ *gzp;
  *g1p= *g3* zz;
  
  if ( ira != 1)
  {
    *g3=( *g3+ *gzp)* rh;
    *gzp = - zz* c1* gz;
  
    if ( rh <= 1.0e-10)
    {
      *g2=0.;
      *g2p=0.;
      return;
    }
  
    *g2= *g2/ rh;
    *g2p= *g2p* zz/ rh;
    return;
  } /* if ( ira != 1) */
  
  t2=.5* a;
  *g2 = - t2* c1* gz;
  *g2p= t2* gz* c2/ r2;
  *g3= rh2* *g2p- a* gz* c1;
  *g2p= *g2p* zz;
  *gzp = - zz* c1* gz;
}

/*-----------------------------------------------------------------------*/

/*
  hfk computes the H-field of a uniform current
  filament by numerical integration
*/
void nec_context::hfk( nec_float el1, nec_float el2, nec_float rhk,
    nec_float zpkx, nec_float *sgr, nec_float *sgi )
{
  int nx = 1, nma = 65536, nts = 4;
  int ns, nt;
  bool flag = true;
  nec_float rx = 1.0e-4;
  nec_float z, ze, s, ep, zend, dz=0., zp, dzot=0., t00r, g1r, g5r, t00i;
  nec_float g1i, g5i, t01r, g3r, t01i, g3i, t10r, t10i, te1i, te1r, t02r;
  nec_float g2r, g4r, t02i, g2i, g4i, t11r, t11i, t20r, t20i, te2i, te2r;
  
  zpka= zpkx;
  rhks= rhk* rhk;
  z= el1;
  ze= el2;
  s= ze- z;
  ep= s/(10.* nma);
  zend= ze- ep;
  *sgr=0.0;
  *sgi=0.0;
  ns= nx;
  nt=0;
  gh( z, &g1r, &g1i);
  
  while( true )
  {
    if ( flag )
    {
      dz= s/ ns;
      zp= z+ dz;
      
      if ( zp > ze )
      {
        dz= ze- z;
        if ( fabs(dz) <= ep )
        {
          *sgr= *sgr* rhk*.5;
          *sgi= *sgi* rhk*.5;
          return;
        }
      }
      
      dzot= dz*.5;
      zp= z+ dzot;
      gh( zp, &g3r, &g3i);
      zp= z+ dz;
      gh( zp, &g5r, &g5i);
    } /* if ( flag ) */
    
    t00r=( g1r+ g5r)* dzot;
    t00i=( g1i+ g5i)* dzot;
    t01r=( t00r+ dz* g3r)*0.5;
    t01i=( t00i+ dz* g3i)*0.5;
    t10r=(4.0* t01r- t00r)/3.0;
    t10i=(4.0* t01i- t00i)/3.0;
    
    test( t01r, t10r, &te1r, t01i, t10i, &te1i, 0.);
    if ( (te1i <= rx) && (te1r <= rx) )
    {
      *sgr= *sgr+ t10r;
      *sgi= *sgi+ t10i;
      nt += 2;
      
      z += dz;
      if ( z >= zend)
      {
        *sgr= *sgr* rhk*.5;
        *sgi= *sgi* rhk*.5;
        return;
      }
      
      g1r= g5r;
      g1i= g5i;
      if ( nt >= nts)
      {
        if ( ns > nx)
        {
          ns= ns/2;
          nt=1;
        }
      }
      flag = true;
      continue;
    } /* if ( (te1i <= rx) && (te1r <= rx) ) */
    
    zp= z+ dz*0.25;
    gh( zp, &g2r, &g2i);
    zp= z+ dz*0.75;
    gh( zp, &g4r, &g4i);
    t02r=( t01r+ dzot*( g2r+ g4r))*0.5;
    t02i=( t01i+ dzot*( g2i+ g4i))*0.5;
    t11r=(4.0* t02r- t01r)/3.0;
    t11i=(4.0* t02i- t01i)/3.0;
    t20r=(16.0* t11r- t10r)/15.0;
    t20i=(16.0* t11i- t10i)/15.0;
    
    test( t11r, t20r, &te2r, t11i, t20i, &te2i, 0.);
    if ( (te2i > rx) || (te2r > rx) )
    {
      nt=0;
      if ( ns >= nma)
        m_output.nec_printf( "\n  STEP SIZE LIMITED AT Z= %10.5f", z );
      else
      {
        ns= ns*2;
        dz= s/ ns;
        dzot= dz*0.5;
        g5r= g3r;
        g5i= g3i;
        g3r= g2r;
        g3i= g2i;
      
        flag = false;
        continue;
      }
    } /* if ( (te2i > rx) || (te2r > rx) ) */
    
    *sgr= *sgr+ t20r;
    *sgi= *sgi+ t20i;
    nt++;
    
    z += dz;
    if ( z >= zend)
    {
      *sgr= *sgr* rhk*.5;
      *sgi= *sgi* rhk*.5;
      return;
    }
    
    g1r= g5r;
    g1i= g5i;
    if ( nt >= nts)
      if ( ns > nx)
      {
        ns= ns/2;
        nt=1;
      }
    flag = true;
    
  } /* while( true ) */
}

/*-----------------------------------------------------------------------*/

/* hintg computes the h field of a patch current */
void nec_context::hintg( nec_float xi, nec_float yi, nec_float zi )
{
  nec_float rx, ry, rfl, xymag, pxx, pyy, cth;
  nec_float rz, rsq, r, rk, cr, sr, t1zr, t2zr;
  nec_complex  gam, f1x, f1y, f1z, f2x, f2y, f2z, rrv, rrh;
  
  rx= xi- xj;
  ry= yi- yj;
  rfl = -1.;
  exk=cplx_00();
  eyk=cplx_00();
  ezk=cplx_00();
  exs=cplx_00();
  eys=cplx_00();
  ezs=cplx_00();
  
  for (int ground_loop = 1; ground_loop <= ground.ksymp; ground_loop++ )
  {
    rfl = - rfl;
    rz = zi - zj*rfl;
    rsq = rx*rx + ry*ry + rz*rz;
  
    if ( rsq < 1.0e-20)
      continue;
  
    r = sqrt( rsq );
    rk= two_pi() * r;
    cr= cos( rk);
    sr= sin( rk);
    gam = -( nec_complex(cr,-sr)+rk*nec_complex(sr,cr) )/( four_pi()*rsq*r )* m_s;
    exc= gam* rx;
    eyc= gam* ry;
    ezc= gam* rz;
    t1zr= t1zj* rfl;
    t2zr= t2zj* rfl;
    f1x= eyc* t1zr- ezc* t1yj;
    f1y= ezc* t1xj- exc* t1zr;
    f1z= exc* t1yj- eyc* t1xj;
    f2x= eyc* t2zr- ezc* t2yj;
    f2y= ezc* t2xj- exc* t2zr;
    f2z= exc* t2yj- eyc* t2xj;
  
    if ( ground_loop != 1)
    {
      if (  ground.type_perfect() ) //ground.iperf == 1)
      {
        f1x = - f1x;
        f1y = - f1y;
        f1z = - f1z;
        f2x = - f2x;
        f2y = - f2y;
        f2z = - f2z;
      }
      else
      {
        xymag= sqrt(rx*rx + ry*ry);
        if ( xymag <= 1.0e-6)
        {
          pxx=0.;
          pyy=0.;
          cth=1.;
          rrv=cplx_10();
        }
        else
        {
          pxx = - ry/ xymag;
          pyy= rx/ xymag;
          cth= rz/ r;
          rrv= sqrt(1.0 - ground.get_zrati_sqr() *(1.0 - cth* cth));
        } /* if ( xymag <= 1.0e-6) */
      
        rrh= ground.zrati* cth;
        rrh=( rrh- rrv)/( rrh+ rrv);
        rrv= ground.zrati* rrv;
        rrv = -( cth- rrv)/( cth+ rrv);
        
        gam=( f1x* pxx+ f1y* pyy)*( rrv- rrh);
        f1x= f1x* rrh+ gam* pxx;
        f1y= f1y* rrh+ gam* pyy;
        f1z= f1z* rrh;
        gam=( f2x* pxx+ f2y* pyy)*( rrv- rrh);
        f2x= f2x* rrh+ gam* pxx;
        f2y= f2y* rrh+ gam* pyy;
        f2z= f2z* rrh;
      } /* if (  ground.iperf == 1) */
    } /* if ( ground_loop != 1) */
  
    exk += f1x;
    eyk += f1y;
    ezk += f1z;
    exs += f2x;
    eys += f2y;
    ezs += f2z;
  
  } /* for( ground_loop = 1; ground_loop <= ground.ksymp; ground_loop++ ) */
}

/*-----------------------------------------------------------------------*/

/*!\brief hsfld computes the h field for constant, sine, and cosine current on a segment including ground effects.
*/
void nec_context::hsfld( nec_float xi, nec_float yi, nec_float zi, nec_float ai )
{
  nec_float xij, yij, rfl, salpr, zij, zp, rhox, rhoy, rhoz, rh, phx;
  nec_float phy, phz, rmag, xymag, xspec, yspec, rhospc, px, py, cth;
  nec_complex hpk, hps, hpc, qx, qy, qz, rrv, rrh, zratx;
  
  xij= xi- xj;
  yij= yi- yj;
  rfl = -1.;
  
  for (int ground_loop = 0; ground_loop < ground.ksymp; ground_loop++ )
  {
    rfl = - rfl;
    salpr= salpj* rfl;
    zij= zi- rfl* zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;
    rh= sqrt( rhox* rhox+ rhoy* rhoy+ rhoz* rhoz+ ai* ai);
  
    if ( rh <= 1.0e-10)
    {
      exk=0.;
      eyk=0.;
      ezk=0.;
      exs=0.;
      eys=0.;
      ezs=0.;
      exc=0.;
      eyc=0.;
      ezc=0.;
      continue;
    }
  
    rhox= rhox/ rh;
    rhoy= rhoy/ rh;
    rhoz= rhoz/ rh;
    phx= sabj* rhoz- salpr* rhoy;
    phy= salpr* rhox- cabj* rhoz;
    phz= cabj* rhoy- sabj* rhox;
  
    hsflx( m_s, rh, zp, &hpk, &hps, &hpc);
  
    if ( ground_loop == 1 )
    {
      if (  false == ground.type_perfect() )
      {
        zratx= ground.zrati;
        rmag= sqrt( zp* zp+ rh* rh);
        xymag= sqrt( xij* xij+ yij* yij);
      
        /* set parameters for radial wire ground screen. */
        if (  ground.get_radial_wire_count() != 0)
        {
          xspec=( xi* zj+ zi* xj)/( zi+ zj);
          yspec=( yi* zj+ zi* yj)/( zi+ zj);
          rhospc= sqrt( xspec*xspec + yspec*yspec + ground.t2*ground.t2);
        
          if ( rhospc <= ground.get_radial_wire_length_wavelengths())
          {
            rrv = ground.m_t1 * rhospc* log(rhospc/ ground.t2);
            zratx = ( rrv* ground.zrati)/( em::impedance() * ground.zrati + rrv);
          }
        }
      
        /* calculation of reflection coefficients when ground is specified. */
        if ( xymag <= 1.0e-6)
        {
          px=0.;
          py=0.;
          cth=1.;
          rrv=cplx_10();
        }
        else
        {
          px = - yij/ xymag;
          py= xij/ xymag;
          cth= zij/ rmag;
          rrv= sqrt(1.- zratx* zratx*(1.- cth* cth));
        }
      
        rrh= zratx* cth;
        rrh = -( rrh- rrv)/( rrh+ rrv);
        rrv= zratx* rrv;
        rrv=( cth- rrv)/( cth+ rrv);
        qy=( phx* px+ phy* py)*( rrv- rrh);
        qx= qy* px+ phx* rrh;
        qy= qy* py+ phy* rrh;
        qz= phz* rrh;
        exk= exk- hpk* qx;
        eyk= eyk- hpk* qy;
        ezk= ezk- hpk* qz;
        exs= exs- hps* qx;
        eys= eys- hps* qy;
        ezs= ezs- hps* qz;
        exc= exc- hpc* qx;
        eyc= eyc- hpc* qy;
        ezc= ezc- hpc* qz;
        continue;
      
      } /* if (  ground.iperf != 1 ) */
      
      exk= exk- hpk* phx;
      eyk= eyk- hpk* phy;
      ezk= ezk- hpk* phz;
      exs= exs- hps* phx;
      eys= eys- hps* phy;
      ezs= ezs- hps* phz;
      exc= exc- hpc* phx;
      eyc= eyc- hpc* phy;
      ezc= ezc- hpc* phz;
      continue;
    
    } /* if ( ground_loop == 1 ) */
  
    exk= hpk* phx;
    eyk= hpk* phy;
    ezk= hpk* phz;
    exs= hps* phx;
    eys= hps* phy;
    ezs= hps* phz;
    exc= hpc* phx;
    eyc= hpc* phy;
    ezc= hpc* phz;
  
  } /* for( ground_loop = 0; ground_loop < ground.ksymp; ground_loop++ ) */
}

/*-----------------------------------------------------------------------*/

/* calculates h field of sine cosine, and constant current of segment */
void nec_context::hsflx( nec_float s, nec_float rh, nec_float zpx,
    nec_complex *hpk, nec_complex *hps,
    nec_complex *hpc )
{
  nec_float r1, r2, zp, z2a, hss, dh, z1;
  nec_float rhz, dk, cdk, sdk, hkr, hki, rh2;
  nec_complex fjk, ekr1, ekr2, t1, t2, cons;
  nec_float pi8 = pi() * 8.0;
  fjk = -two_pi_j();
  if ( rh >= 1.0e-10)
  {
    if ( zpx >= 0.0)
    {
      zp= zpx;
      hss=1.;
    }
    else
    {
      zp = - zpx;
      hss = -1.;
    }
  
    dh= 0.5* s;
    z1= zp+ dh;
    z2a= zp- dh;
    if ( z2a >= 1.0e-7)
      rhz= rh/ z2a;
    else
      rhz=1.;
  
    dk= two_pi() * dh;
    cdk= cos( dk);
    sdk= sin( dk);
    hfk(- dk, dk, rh* two_pi(), zp* two_pi(), &hkr, &hki);
    *hpk= nec_complex( hkr, hki);
  
    if ( rhz >= 1.0e-3)
    {
      rh2= rh* rh;
      r1= sqrt( rh2+ z1* z1);
      r2= sqrt( rh2+ z2a* z2a);
      ekr1= exp( fjk* r1);
      ekr2= exp( fjk* r2);
      t1= z1* ekr1/ r1;
      t2= z2a* ekr2/ r2;
      *hps=( cdk*( ekr2- ekr1)- cplx_01()* sdk*( t2+ t1))* hss;
      *hpc = - sdk*( ekr2+ ekr1)- cplx_01()* cdk*( t2- t1);
      cons = - cplx_01()/(2.0 * two_pi() * rh);
      *hps= cons* *hps;
      *hpc= cons* *hpc;
      return;
    
    } /* if ( rhz >= 1.0e-3) */
  
    ekr1= nec_complex( cdk, sdk)/( z2a* z2a);
    ekr2= nec_complex( cdk,- sdk)/( z1* z1);
    t1= two_pi()*(1.0/z1 - 1.0/z2a);
    t2= exp( fjk* zp)* rh/ pi8;
    *hps= t2*( t1+( ekr1+ ekr2)* sdk)* hss;
    *hpc= t2*(- cplx_01()* t1+( ekr1- ekr2)* cdk);
    return;
  
  } /* if ( rh >= 1.0e-10) */
  
  *hps=cplx_00();
  *hpc=cplx_00();
  *hpk=cplx_00();
}



/*
  intx performs numerical integration of exp(jkr)/r by the method of 
  variable interval width Romberg integration.  The integrand value
  is supplied by subroutine gf.
*/
void nec_context::intx( nec_float el1, nec_float el2, nec_float b, int ij, nec_float *sgr, nec_float *sgi)
{
  int ns, nt;
  int nx = 1, nma = 65536, nts = 4;
  bool flag = true;
  nec_float z, s, ze, fnm, ep, zend, fns, dz=0., zp, dzot=0., t00r, g1r, g5r, t00i;
  nec_float g1i, g5i, t01r, g3r, t01i, g3i, t10r, t10i, te1i, te1r, t02r;
  nec_float g2r, g4r, t02i, g2i, g4i, t11r, t11i, t20r, t20i, te2i, te2r;
  nec_float rx = 1.0e-4;

  z= el1;
  ze= el2;
  if ( ij == 0)
    ze=0.;
  s= ze- z;
  fnm= nma;
  ep= s/(10.* fnm);
  zend= ze- ep;
  *sgr=0.;
  *sgi=0.;
  ns= nx;
  nt=0;
  gf( z, &g1r, &g1i);

  while( true )
  {
    if ( flag )
    {
      fns= ns;
      dz= s/ fns;
      zp= z+ dz;

      if ( zp > ze)
      {
  dz= ze- z;
  if ( fabs(dz) <= ep)
  {
    /* add contribution of near singularity for diagonal term */
    if (ij == 0)
    {
      *sgr=2.*( *sgr+ log(( sqrt( b* b+ s* s)+ s)/ b));
      *sgi=2.* *sgi;
    }
    return;
  }

      } /* if ( zp > ze) */

      dzot= dz*.5;
      zp= z+ dzot;
      gf( zp, &g3r, &g3i);
      zp= z+ dz;
      gf( zp, &g5r, &g5i);

    } /* if ( flag ) */

    t00r=( g1r+ g5r)* dzot;
    t00i=( g1i+ g5i)* dzot;
    t01r=( t00r+ dz* g3r)*0.5;
    t01i=( t00i+ dz* g3i)*0.5;
    t10r=(4.0* t01r- t00r)/3.0;
    t10i=(4.0* t01i- t00i)/3.0;

    /* test convergence of 3 point romberg result. */
    test( t01r, t10r, &te1r, t01i, t10i, &te1i, 0.);
    if ( (te1i <= rx) && (te1r <= rx) )
    {
      *sgr= *sgr+ t10r;
      *sgi= *sgi+ t10i;
      nt += 2;

      z += dz;
      if ( z >= zend)
      {
  /* add contribution of near singularity for diagonal term */
  if (ij == 0)
  {
    *sgr=2.*( *sgr+ log(( sqrt( b* b+ s* s)+ s)/ b));
    *sgi=2.* *sgi;
  }
  return;
      }

      g1r= g5r;
      g1i= g5i;
      if ( nt >= nts)
  if ( ns > nx)
  {
    /* Double step size */
    ns= ns/2;
    nt=1;
  }
      flag = true;
      continue;

    } /* if ( (te1i <= rx) && (te1r <= rx) ) */

    zp= z+ dz*0.25;
    gf( zp, &g2r, &g2i);
    zp= z+ dz*0.75;
    gf( zp, &g4r, &g4i);
    t02r=( t01r+ dzot*( g2r+ g4r))*0.5;
    t02i=( t01i+ dzot*( g2i+ g4i))*0.5;
    t11r=(4.0* t02r- t01r)/3.0;
    t11i=(4.0* t02i- t01i)/3.0;
    t20r=(16.0* t11r- t10r)/15.0;
    t20i=(16.0* t11i- t10i)/15.0;

    /* test convergence of 5 point romberg result. */
    test( t11r, t20r, &te2r, t11i, t20i, &te2i, 0.);
    if ( (te2i > rx) || (te2r > rx) )
    {
      nt=0;
      if ( ns >= nma)
  m_output.nec_printf( "\n  STEP SIZE LIMITED AT Z= %10.5f", z );
      else
      {
  /* halve step size */
  ns= ns*2;
  fns= ns;
  dz= s/ fns;
  dzot= dz*0.5;
  g5r= g3r;
  g5i= g3i;
  g3r= g2r;
  g3i= g2i;

  flag = false;
  continue;
      }

    } /* if ( (te2i > rx) || (te2r > rx) ) */

    *sgr= *sgr+ t20r;
    *sgi= *sgi+ t20i;
    nt++;

    z += dz;
    if ( z >= zend)
    {
      /* add contribution of near singularity for diagonal term */
      if (ij == 0)
      {
  *sgr=2.*( *sgr+ log(( sqrt( b* b+ s* s)+ s)/ b));
  *sgi=2.* *sgi;
      }
      return;
    }

    g1r= g5r;
    g1i= g5i;
    if ( nt >= nts)
      if ( ns > nx)
      {
  /* Double step size */
  ns= ns/2;
  nt=1;
      }
    flag = true;

  } /* while( true ) */

}



/*! \brief nefld computes the near field at specified points in space after the structure currents have been computed.
  \param xob,yob,zob The field evaluation points.
*/
void nec_context::nefld( nec_float xob, nec_float yob, nec_float zob,
    nec_complex *ex, nec_complex *ey, nec_complex *ez )
{
  int i, ix, ipr, jc;
  nec_float zp, ax;
  nec_complex acx, bcx, ccx;
  
  *ex=cplx_00();
  *ey=cplx_00();
  *ez=cplx_00();
  ax=0.0;
  
  int n = m_geometry->n_segments;
  for( i = 0; i < n; i++ )
  {
    xj= xob- m_geometry->x[i];
    yj= yob- m_geometry->y[i];
    zj= zob- m_geometry->z[i];
    zp= m_geometry->cab[i]* xj+ m_geometry->sab[i]* yj+ m_geometry->salp[i]* zj;
  
    if ( fabs( zp) > 0.5001* m_geometry->segment_length[i])
      continue;
  
    zp= xj* xj+ yj* yj+ zj* zj- zp* zp;
    xj= m_geometry->segment_radius[i];
  
    if ( zp > 0.9* xj* xj)
      continue;
  
    ax= xj;
    break;
  } /* for( i = 0; i < n; i++ ) */

  for( i = 0; i < n; i++ )
  {
    ix = i+1;
    m_s = m_geometry->segment_length[i];
    m_b = m_geometry->segment_radius[i];
    xj= m_geometry->x[i];
    yj= m_geometry->y[i];
    zj= m_geometry->z[i];
    cabj= m_geometry->cab[i];
    sabj= m_geometry->sab[i];
    salpj= m_geometry->salp[i];
  
    if ( m_use_exk == true)
    {
    ipr= m_geometry->icon1[i];
  
    if (ipr > PCHCON) ind1 = 2;
    else if ( ipr < 0 )
    {
      ipr = -ipr;
      int iprx = ipr-1;
    
      if ( -m_geometry->icon1[iprx] != ix )
        ind1=2;
      else
      {
        ind1 = m_geometry->test_ek_approximation(i,iprx);
      }
    } /* if ( ipr < 0 ) */
    else if ( ipr == 0 )
      ind1=1;
    else
    {
      int iprx = ipr-1;
  
      if ( ipr != ix )
      {
        if ( m_geometry->icon2[iprx] != ix )
          ind1=2;
        else
        {
          ind1 = m_geometry->test_ek_approximation(i,iprx);
        }
      } /* if ( ipr != ix ) */
      else
      {
        if ( (cabj* cabj+ sabj* sabj) > 1.0e-8)
          ind1=2;
        else
          ind1=0;
      }
    } /* else */
  
    ipr= m_geometry->icon2[i];
  
    if (ipr > PCHCON) ind2 = 2;
    else if ( ipr < 0 )
    {
      ipr = -ipr;
      int iprx = ipr-1;
    
      if ( -m_geometry->icon2[iprx] != ix )
        ind1=2;
      else
      {
        ind1 = m_geometry->test_ek_approximation(i,iprx);
      }
    } /* if ( ipr < 0 ) */
    else if ( ipr == 0 )
      ind2=1;
    else
    {
      int iprx = ipr-1;
  
      if ( ipr != ix )
      {
        if ( m_geometry->icon1[iprx] != ix )
          ind2=2;
        else
        {
          ind2 = m_geometry->test_ek_approximation(i,iprx);
        }
      } /* if ( ipr != (i+1) ) */
      else
      {
        if ( (cabj* cabj+ sabj* sabj) > 1.0e-8)
          ind1=2;
        else
          ind1=0;
      }
    } /* else */
  
    } /* if ( m_use_exk == true) */
  
    efld( xob, yob, zob, ax, true);
    acx= nec_complex( air[i], aii[i]);
    bcx= nec_complex( bir[i], bii[i]);
    ccx= nec_complex( cir[i], cii[i]);
    *ex += exk* acx+ exs* bcx+ exc* ccx;
    *ey += eyk* acx+ eys* bcx+ eyc* ccx;
    *ez += ezk* acx+ ezs* bcx+ ezc* ccx;
  
  } /* for( i = 0; i < n; i++ ) */
  
  
  jc= n-1;
  for( i = 0; i < m_geometry->m; i++ )
  {
    m_s= m_geometry->pbi[i];
    xj= m_geometry->px[i];
    yj= m_geometry->py[i];
    zj= m_geometry->pz[i];
    t1xj= m_geometry->t1x[i];
    t1yj= m_geometry->t1y[i];
    t1zj= m_geometry->t1z[i];
    t2xj= m_geometry->t2x[i];
    t2yj= m_geometry->t2y[i];
    t2zj= m_geometry->t2z[i];
    jc += 3;
    acx= t1xj* current_vector[jc-2]+ t1yj* current_vector[jc-1]+ t1zj* current_vector[jc];
    bcx= t2xj* current_vector[jc-2]+ t2yj* current_vector[jc-1]+ t2zj* current_vector[jc];
  
    {
      unere( xob, yob, zob, false);
      *ex += acx* exk+ bcx* exs;
      *ey += acx* eyk+ bcx* eys;
      *ez += acx* ezk+ bcx* ezs;
    }
    if (ground.present())
    {
      unere( xob, yob, zob, true);
      *ex += acx* exk+ bcx* exs;
      *ey += acx* eyk+ bcx* eys;
      *ez += acx* ezk+ bcx* ezs;
    }
/*    for( ipa = 0; ipa < ground.ksymp; ipa++ )
    {
      int ipgnd= ipa+1;
      unere( xob, yob, zob, ipgnd == 2);
      *ex= *ex+ acx* exk+ bcx* exs;
      *ey= *ey+ acx* eyk+ bcx* eys;
      *ez= *ez+ acx* ezk+ bcx* ezs;
    } */
  } /* for( i = 0; i < m; i++ ) */
}
  
/*-----------------------------------------------------------------------*/

#include <sstream>
/* subroutine netwk solves for structure currents for a given */
/* excitation including the effect of non-radiating networks if */
/* present. */
void nec_context::netwk( complex_array& in_cm, int_array& in_ip,
    complex_array& einc )
  /*
    The parameters cmb, cmc, cmd appear to be Numerical Green's Functions 
    Related. N.G.F are not implemented in this version of NEC.
    
    The evidence for the above opinion is from the NEC-PC code that does
    implement the N.G.F.
  */
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
    rhs.resize( m_geometry->n_plus_3m ); // this should probably be ndimn!
  
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
          asmx= m_geometry->segment_length[isc1];
        
          vector_fill(rhs,0,neqt,cplx_00());
        
          rhs[isc1] = cplx_10();
          solves( in_cm, in_ip, rhs, neq, 1, m_geometry->np, m_geometry->n_segments, m_geometry->mp, m_geometry->m, nop, symmetry_array);
          m_geometry->get_current_coefficients(_wavelength, rhs, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
        
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
        m_output.nec_printf( "\n\n"
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
      cmn.setConstant(cplx_00());
      rhnx.setConstant(cplx_00());
      
      nteq=0;
      ntsc=0;
  
      /*  sort network and source data and
        assign equation numbers to segments */
      for (int j = 0; j < network_count; j++ )
      {
        nseg1= iseg1[j];
        nseg2= iseg2[j];
      
        // Calculate the admittance yXX for this network element
        if ( ntyp[j] <= 1) // if not transmission line
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
          y22r= two_pi() * x11i[j]/ _wavelength;
          y12r=0.;
          y12i=1./( x11r[j]* sin( y22r));
          y11r= x12r[j];
          y11i = - y12i* cos( y22r);
          y22r= x22r[j];
          y22i= y11i+ x22i[j];
          y11i= y11i+ x12i[j];
        
          if ( ntyp[j] != 2) // ntype == 3, must be crossed element
          {
            y12r = -y12r;
            y12i = -y12i;
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
          isc1 = -1;
      
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
          isc2 = -1;
      
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
    
        /*
          Fill network equation matrix and right hand side vector with
          network short-circuit admittance matrix coefficients.
        */
        if ( isc1 == -1)
        {
          cmn[irow1+irow1*ndimn] -= nec_complex( y11r, y11i)* m_geometry->segment_length[nseg1-1];
          cmn[irow1+irow2*ndimn] -= nec_complex( y12r, y12i)* m_geometry->segment_length[nseg1-1];
        }
        else
        {
          rhnx[irow1] += nec_complex( y11r, y11i)* source_voltage_array[isc1]/_wavelength;
          rhnx[irow2] += nec_complex( y12r, y12i)* source_voltage_array[isc1]/_wavelength;
        }
      
        if ( isc2 == -1)
        {
          cmn[irow2+irow2*ndimn] -= nec_complex( y22r, y22i)* m_geometry->segment_length[nseg2-1];
          cmn[irow2+irow1*ndimn] -= nec_complex( y12r, y12i)* m_geometry->segment_length[nseg2-1];
        }
        else
        {
          rhnx[irow1] += nec_complex( y12r, y12i)* source_voltage_array[isc2]/_wavelength;
          rhnx[irow2] += nec_complex( y22r, y22i)* source_voltage_array[isc2]/_wavelength;
        }
      } /* for( j = 0; j < network_count; j++ ) */
  
      /*  add interaction matrix admittance
        elements to network equation matrix */
      for( i = 0; i < nteq; i++ )
      {
        vector_fill(rhs, 0, neqt, cplx_00());
        
        irow1= nteqa[i]-1;
        rhs[irow1]=cplx_10();
        solves( in_cm, in_ip, rhs, neq, 1, m_geometry->np, m_geometry->n_segments, m_geometry->mp, m_geometry->m, nop, symmetry_array);
        m_geometry->get_current_coefficients(_wavelength, rhs, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
        
        for (int j = 0; j < nteq; j++ )
        {
          irow1= nteqa[j]-1;
          cmn[i+j*ndimn] += rhs[irow1];
        }
      } /* for( i = 0; i < nteq; i++ ) */
    
      /* factor network equation matrix */
      lu_decompose(m_output, nteq, cmn, ipnt, ndimn);
      
    } /* if ( network_count != 0) */
  } /* if ( ntsol != 0) */

  if (0 == network_count)
  {
    /* solve for currents when no networks are present */
    solves( in_cm, in_ip, einc, neq, 1, m_geometry->np, m_geometry->n_segments, m_geometry->mp, m_geometry->m, nop, symmetry_array);
    m_geometry->get_current_coefficients(_wavelength, einc, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
    ntsc=0;
  }
  else // if ( network_count != 0)
  {
    /*
      Add to network equation right hand side
      the terms due to element interactions
    */
    for (i = 0; i < neqt; i++)
      rhs[i]= einc[i];
  
    solves( in_cm, in_ip, rhs, neq, 1, m_geometry->np, m_geometry->n_segments, m_geometry->mp, m_geometry->m, nop, symmetry_array);
    m_geometry->get_current_coefficients(_wavelength, rhs, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);
  
    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      rhnt[i]= rhnx[i]+ rhs[irow1];
    }

    /* solve network equations */
    solve( nteq, cmn, ipnt, rhnt, ndimn);
  
    /*
      Add fields due to network voltages to electric fields
      applied to structure and solve for induced current
    */
    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      einc[irow1] -= rhnt[i];
    }
  
    solves( in_cm, in_ip, einc, neq, 1, m_geometry->np, m_geometry->n_segments, m_geometry->mp, m_geometry->m, nop, symmetry_array);
    m_geometry->get_current_coefficients(_wavelength, einc, air, aii, bir, bii, cir, cii, vqds, nqds, iqds);


    nec_structure_excitation* seo = new nec_structure_excitation();
    
    for( i = 0; i < nteq; i++ )
    {
      int segment_number = nteqa[i];
      int segment_index = segment_number-1;
      nec_complex voltage = rhnt[i]* m_geometry->segment_length[segment_index]* _wavelength;
      nec_complex current = einc[segment_index]* _wavelength;
      nec_float power = em::power(voltage,current);
      network_power_loss = network_power_loss - power;
      
      int segment_tag = m_geometry->segment_tags[segment_index];
      seo->add(segment_number, segment_tag, voltage, current, power);
    }
    
    for ( i = 0; i < ntsc; i++ )
    {
      int segment_number = ntsca[i];
      int segment_index = segment_number-1;
      nec_complex voltage = vsrc[i];
      nec_complex current = einc[segment_index]* _wavelength;
      nec_float power = em::power(voltage,current);
      network_power_loss = network_power_loss - power;
      
      int segment_tag = m_geometry->segment_tags[segment_index];
      seo->add(segment_number, segment_tag, voltage, current, power);
    } /* for( i = 0; i < ntsc; i++ ) */
    
    if ( nprint == 0)
    {
      std::stringstream ss;
      seo->write_to_file(ss);
      m_output.line(ss.str().c_str());
    }
    
    seo->set_frequency(freq_mhz/(1.e-6));
    m_results.add(seo);

      
/*    if ( nprint == 0)
    {
      m_output.nec_printf( "\n\n\n"
        "                          "
        "--------- STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS --------" );
    
      m_output.nec_printf( "\n"
        "  TAG   SEG       VOLTAGE (VOLTS)          CURRENT (AMPS)        "
        " IMPEDANCE (OHMS)       ADMITTANCE (MHOS)     POWER\n"
        "  No:   No:     REAL      IMAGINARY     REAL      IMAGINARY    "
        " REAL      IMAGINARY     REAL      IMAGINARY   (WATTS)" );
    }

    for( i = 0; i < nteq; i++ )
    {
      int segment_number = nteqa[i];
      int segment_index = segment_number-1;
      nec_complex voltage = rhnt[i]* m_geometry->segment_length[segment_index]* _wavelength;
      nec_complex current = einc[segment_index]* _wavelength;
      nec_complex admittance = current / voltage;
      nec_complex impedance = voltage / current;
      int segment_tag = m_geometry->segment_tags[segment_number];
      nec_float power = em::power(voltage,current);
      network_power_loss = network_power_loss - power;
      
      if ( nprint == 0)
        m_output.nec_printf( "\n"
          " %4d %5d %11.4E %11.4E %11.4E %11.4E"
          " %11.4E %11.4E %11.4E %11.4E %11.4E",
          segment_tag, segment_number, real(voltage), imag(voltage), real(current), imag(current),
          real(impedance), imag(impedance), real(admittance), imag(admittance), power );
    }

    for( i = 0; i < ntsc; i++ )
    {
      irow1= ntsca[i]-1;
      vlt= vsrc[i];
      cux= einc[irow1]* _wavelength;
      ymit= cux/ vlt;
      zped= vlt/ cux;
      irow2= m_geometry->segment_tags[irow1];
      
      nec_float pwr= em::power(vlt,cux);
      network_power_loss= network_power_loss- pwr;
    
      if ( nprint == 0)
        m_output.nec_printf( "\n"
          " %4d %5d %11.4E %11.4E %11.4E %11.4E"
          " %11.4E %11.4E %11.4E %11.4E %11.4E",
          irow2, irow1+1, real(vlt), imag(vlt), real(cux), imag(cux),
          real(zped), imag(zped), real(ymit), imag(ymit), pwr );
    } // for( i = 0; i < ntsc; i++ )
*/
  } // if ( network_count != 0)

  if ( (voltage_source_count+nvqd) == 0)
    return;
  
  
  // Create an antenna_input results object to hold the antenna input results.
  nec_antenna_input* antenna_input = new nec_antenna_input();
  antenna_input->set_frequency(freq_mhz/(1.e-6));
  m_results.add(antenna_input);
  DEBUG_TRACE("Creating new antenna input")
/*  m_output.end_section();
  m_output.nec_printf( 
    "                        "
    "--------- ANTENNA INPUT PARAMETERS ---------" );
  
  m_output.nec_printf( "\n"
    "  TAG   SEG       VOLTAGE (VOLTS)         "
    "CURRENT (AMPS)         IMPEDANCE (OHMS)    "
    "    ADMITTANCE (MHOS)     POWER\n"
    "  NO.   NO.     REAL      IMAGINARY"
    "     REAL      IMAGINARY     REAL      "
    "IMAGINARY    REAL       IMAGINARY   (WATTS)" );
*/  
  for( i = 0; i < voltage_source_count; i++ )
  {
    int segment_index = source_segment_array[i]-1;
    nec_complex voltage = source_voltage_array[i];
    nec_complex current = einc[segment_index] * _wavelength;
    
    bool add_as_network_loss = false;
    
    if( ntsc != 0)
    {
      int j = 0;
      for (j = 0; j < ntsc; j++ )
        if (ntsca[j] == segment_index+1)
          break;
      
      int row_index = ndimn - (j+1);
      int row_offset = row_index*ndimn;
      
      nec_complex temp = rhnx[row_index];
      for (int k = 0; k < nteq; k++ )
        temp -= cmn[k+row_offset]*rhnt[k];
        
      current = (einc[segment_index] + temp) * _wavelength;
      add_as_network_loss = true;
    }
    
/*
  The following Code provided by Neoklis Kyriazis as a replacement for
  the broken code below.

  if( ntsc == 0)
  {
    cux= einc[isc1]* wlam;
    irow1=0;
  }
  else
  {
    for( j = 0; j < ntsc; j++ )
      if( ntsca[j] == isc1+1)
        break;

    irow1= ndimn- (j+1);
    cux= rhnx[irow1];
    for( j = 0; j < nteq; j++ )
      cux -= cmn[j+irow1*ndimn]*rhnt[j];
    cux=(einc[isc1]+ cux)* wlam;
    irow1++;

  } // if( ntsc == 0)
    
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
          
        current = (temp + einc[segment_index])* _wavelength;
        add_as_network_loss = true;
          
#warning "This loop is messed up. The j is inside another j loop"
        // I have removed the j from the "for (int k = 0; k < nteq; k++ )" loop 
        // and placed this"j=nteq" statement here.
        j = nteq;
      }
    }
*/
    nec_complex admittance = current / voltage;
    nec_complex impedance = voltage / current;
    nec_float power = em::power(voltage,current);
    
    if ( add_as_network_loss )
      network_power_loss += power;
      
    input_power += power;
    
    int segment_tag = m_geometry->segment_tags[segment_index];
  
    antenna_input->set_input(
      segment_tag, segment_index+1,
      voltage, current, impedance, admittance, power);
    
/*    m_output.nec_printf(  "\n"
      " %4d %5d %11.4E %11.4E %11.4E %11.4E"
      " %11.4E %11.4E %11.4E %11.4E %11.4E",
      segment_tag, segment_index+1, real(voltage), imag(voltage), real(current), imag(current),
      real(impedance), imag(impedance), real(admittance), imag(admittance), power ); */
    
  } /* for( i = 0; i < voltage_source_count; i++ ) */

  
  for( i = 0; i < nvqd; i++ )
  {
/*
    isc1= ivqd[i]-1;
    vlt= vqd[i];
    cux= cmplx( air[isc1], aii[isc1]);
    ymit= cmplx( bir[isc1], bii[isc1]);
    zped= cmplx( cir[isc1], cii[isc1]);
    pwr= si[isc1]* TP*.5;
    cux=( cux- ymit* sin( pwr)+ zped* cos( pwr))* wlam;
    ymit= cux/ vlt;
    zped= vlt/ cux;
    pwr=.5* creal( vlt* conj( cux));
    pin= pin+ pwr;
    irow2= itag[isc1];

    fprintf( output_fp,  "\n"
  " %4d %5d %11.4E %11.4E %11.4E %11.4E"
  " %11.4E %11.4E %11.4E %11.4E %11.4E",
  irow2, isc1+1, creal(vlt), cimag(vlt), creal(cux), cimag(cux),
  creal(zped), cimag(zped), creal(ymit), cimag(ymit), pwr );
*/
    int segment_index = ivqd[i]-1;
    nec_complex voltage = vqd[i];
    
    nec_complex _ai( air[segment_index], aii[segment_index]);
    nec_complex _bi( bir[segment_index], bii[segment_index]);
    nec_complex _ci( cir[segment_index], cii[segment_index]);
    
    // segment length is measured in wavelengths. The phase is therefore the length in wavelengths
    // multiplied by pi().// TCAM CHANGED TO pi() (from TP*.5)!!
    nec_float segment_length_phase = m_geometry->segment_length[segment_index] * pi(); 
    
    nec_complex current = (_ai - _bi* sin(segment_length_phase) + _ci * cos(segment_length_phase)) * _wavelength;
    
    nec_complex admittance = current / voltage;
    nec_complex impedance = voltage / current;
    nec_float power = em::power(voltage,current);
    
    input_power += power;
    
    int segment_tag = m_geometry->segment_tags[segment_index];
  
    antenna_input->set_input(
      segment_tag, segment_index+1,
      voltage, current, impedance, admittance, power);
    
/*    m_output.nec_printf(  "\n"
      " %4d %5d %11.4E %11.4E %11.4E %11.4E"
      " %11.4E %11.4E %11.4E %11.4E %11.4E",
      segment_tag, segment_index+1, real(voltage), imag(voltage), real(current), imag(current),
      real(impedance), imag(impedance), real(admittance), imag(admittance), power ); */
    
  } /* for( i = 0; i < nvqd; i++ ) */
  
  std::stringstream ss;
  antenna_input->write_to_file(ss);
  m_output.line(ss.str().c_str());
}


/*-----------------------------------------------------------------------*/

/* compute near e or h fields over a range of points */
void nec_context::nfpat( void )
{
  int i, j, kk;
  nec_float znrt, cth=0., sth=0., ynrt, cph=0., sph=0., xnrt, xob, yob;
  nec_float zob, tmp1, tmp2, tmp3;
  /* nec_float tmp4, tmp5, tmp6; */
  nec_complex ex, ey, ez;

/* The printing of near fields is now managed by a dedicated nec_base_result : nec_near_field_pattern*/

/* if ( nfeh != 1)
  {
    m_output.nec_printf(  "\n\n\n"
  "                             "
  "-------- NEAR ELECTRIC FIELDS --------\n"
  "     ------- LOCATION -------     ------- EX ------    ------- EY ------    ------- EZ ------\n"
  "      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n"
  "    METERS    METERS    METERS     VOLTS/M  DEGREES    VOLTS/M   DEGREES     VOLTS/M  DEGREES" );
  }
  else
  {
    m_output.nec_printf(  "\n\n\n"
  "                                   "
  "-------- NEAR MAGNETIC FIELDS ---------\n\n"
  "     ------- LOCATION -------     ------- HX ------    ------- HY ------    ------- HZ ------\n"
  "      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n"
  "    METERS    METERS    METERS      AMPS/M  DEGREES      AMPS/M  DEGREES      AMPS/M  DEGREES" );
  }*/
  
    
  nec_near_field_pattern* nfp = new nec_near_field_pattern(nfeh);  
  nfp->set_frequency(freq_mhz/(1.e-6));
  m_results.add(nfp);
  
  
  znrt= znr- dznr;
        
  for( i = 0; i < nrz; i++ )
  {
    znrt += dznr;
    if ( m_near != 0)
    {
      cth= cos( degrees_to_rad(znrt) );
      sth= sin( degrees_to_rad(znrt) );
    }

    ynrt= ynr- dynr;
    for( j = 0; j < nry; j++ )
    {
      ynrt += dynr;
      if ( m_near != 0)
      {
  cph= cos( degrees_to_rad(ynrt) );
  sph= sin( degrees_to_rad(ynrt) );
      }

      xnrt= xnr- dxnr;
      for( kk = 0; kk < nrx; kk++ )
      {
  xnrt += dxnr;
  if ( m_near != 0)
  {
    xob= xnrt* sth* cph;
    yob= xnrt* sth* sph;
    zob= xnrt* cth;
  }
  else
  {
    xob= xnrt;
    yob= ynrt;
    zob= znrt;
  }

  tmp1= xob/ _wavelength;
  tmp2= yob/ _wavelength;
  tmp3= zob/ _wavelength;

  if ( nfeh != 1)
    nefld( tmp1, tmp2, tmp3, &ex, &ey, &ez);
  else
    nhfld( tmp1, tmp2, tmp3, &ex, &ey, &ez);

  /*no more used */
  
  /* 
  tmp1= abs( ex);
  tmp2= arg_degrees( ex);
  tmp3= abs( ey);
  tmp4= arg_degrees( ey);
  tmp5= abs( ez);
  tmp6= arg_degrees( ez);
  
  m_output.nec_printf( "\n"
      " %9.4f %9.4f %9.4f  %11.4E %7.2f  %11.4E %7.2f  %11.4E %7.2f",
      xob, yob, zob, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );
  */
  
  nfp->set_input(xob, yob, zob, ex, ey, ez);
      
  plot_card.plot_fields(ex, ey, ez, xob, yob, zob);

  
      } /* for( kk = 0; kk < nrx; kk++ ) */

    } /* for( j = 0; j < nry; j++ ) */

  } /* for( i = 0; i < nrz; i++ ) */
  

// write the results to the ordinary NEC output file.
        
std::stringstream ss;
nfp->write_to_file(ss);
m_output.line(ss.str().c_str());  
  
}

/*-----------------------------------------------------------------------*/

/*!\brief Computes the near field at specified points in space after the structure currents have been computed.
*/
void nec_context::nhfld( nec_float xob, nec_float yob, nec_float zob,
    nec_complex *hx, nec_complex *hy, nec_complex *hz )
{
  *hx=cplx_00();
  *hy=cplx_00();
  *hz=cplx_00();
  nec_float ax = 0.0;
  
  int n = m_geometry->n_segments;
  for (int i = 0; i < n; i++ )
  {
    xj= xob- m_geometry->x[i];
    yj= yob- m_geometry->y[i];
    zj= zob- m_geometry->z[i];
    nec_float zp= m_geometry->cab[i]* xj+ m_geometry->sab[i]* yj+ m_geometry->salp[i]* zj;
    
    if ( fabs( zp) > 0.5001* m_geometry->segment_length[i])
      continue;
    
    zp= xj* xj+ yj* yj+ zj* zj- zp* zp;
    xj= m_geometry->segment_radius[i];
    
    if ( zp > 0.9* xj* xj)
      continue;
    
    ax = xj;
    break;
  }
  
  for (int i = 0; i < n; i++ )
  {
    m_s = m_geometry->segment_length[i];
    m_b = m_geometry->segment_radius[i];
    xj= m_geometry->x[i];
    yj= m_geometry->y[i];
    zj= m_geometry->z[i];
    cabj= m_geometry->cab[i];
    sabj= m_geometry->sab[i];
    salpj= m_geometry->salp[i];
    hsfld( xob, yob, zob, ax);
    nec_complex acx= nec_complex( air[i], aii[i]);
    nec_complex bcx= nec_complex( bir[i], bii[i]);
    nec_complex ccx= nec_complex( cir[i], cii[i]);
    *hx += exk* acx+ exs* bcx+ exc* ccx;
    *hy += eyk* acx+ eys* bcx+ eyc* ccx;
    *hz += ezk* acx+ ezs* bcx+ ezc* ccx;
  }
  
  if ( m_geometry->m == 0)
    return;
  
  int jc= n-1;
  for (int i = 0; i < m_geometry->m; i++ )
  {
    m_s = m_geometry->pbi[i];
    xj= m_geometry->px[i];
    yj= m_geometry->py[i];
    zj= m_geometry->pz[i];
    t1xj= m_geometry->t1x[i];
    t1yj= m_geometry->t1y[i];
    t1zj= m_geometry->t1z[i];
    t2xj= m_geometry->t2x[i];
    t2yj= m_geometry->t2y[i];
    t2zj= m_geometry->t2z[i];
    hintg( xob, yob, zob);
    jc += 3;
    nec_complex acx= t1xj* current_vector[jc-2] + t1yj* current_vector[jc-1] + t1zj* current_vector[jc];
    nec_complex bcx= t2xj* current_vector[jc-2] + t2yj* current_vector[jc-1] + t2zj* current_vector[jc];
    *hx= *hx+ acx* exk+ bcx* exs;
    *hy= *hy+ acx* eyk+ bcx* eys;
    *hz= *hz+ acx* ezk+ bcx* ezs;
  }
}
  

/*-----------------------------------------------------------------------*/

/**\brief Integrate over patches at wire connection point 
*/
void nec_context::pcint( nec_float xi, nec_float yi, nec_float zi, nec_float cabi,
    nec_float sabi, nec_float salpi, complex_array& e)
{
  nec_float d, ds, da, gcon, fcon, xxj, xyj, xzj, xs, s1;
  nec_float xss, yss, zss, s2x;
  
  int nint = 10;
  d = sqrt(m_s)*0.5;
  ds=4.* d/ (nec_float) nint;
  da= ds* ds;
  gcon=1./ m_s;
  fcon=1./(2.* two_pi() * d);
  xxj= xj;
  xyj= yj;
  xzj= zj;
  xs = m_s;
  m_s = da;
  s1= d+ ds*.5;
  xss= xj+ s1*( t1xj+ t2xj);
  yss= yj+ s1*( t1yj+ t2yj);
  zss= zj+ s1*( t1zj+ t2zj);
  s1= s1+ d;
  s2x= s1;
  
  nec_complex e1=cplx_00();
  nec_complex e2=cplx_00();
  nec_complex e3=cplx_00();
  nec_complex e4=cplx_00();
  nec_complex e5=cplx_00();
  nec_complex e6=cplx_00();
  nec_complex e7=cplx_00();
  nec_complex e8=cplx_00();
  nec_complex e9=cplx_00();
  
  for (int i1 = 0; i1 < nint; i1++ )
  {
    s1= s1- ds;
    nec_float s2 = s2x;
    xss= xss- ds* t1xj;
    yss= yss- ds* t1yj;
    zss= zss- ds* t1zj;
    xj= xss;
    yj= yss;
    zj= zss;
  
    for(int i2 = 0; i2 < nint; i2++ )
    {
      s2= s2- ds;
      xj= xj- ds* t2xj;
      yj= yj- ds* t2yj;
      zj= zj- ds* t2zj;
      unere( xi, yi, zi, false);
      exk= exk* cabi+ eyk* sabi+ ezk* salpi;
      exs= exs* cabi+ eys* sabi+ ezs* salpi;
      nec_float g1 = (d + s1)*(d + s2)*gcon;
      nec_float g2 = (d - s1)*(d + s2)*gcon;
      nec_float g3 = (d - s1)*(d - s2)*gcon;
      nec_float g4 = (d + s1)*(d - s2)*gcon;
      nec_float f2 = (s1*s1 + s2*s2) * two_pi();
      nec_float f1= s1/f2 - (g1 - g2 - g3 + g4)*fcon;
      f2 = s2/f2 - (g1 + g2 - g3 - g4)*fcon;
      e1 += exk* g1;
      e2 += exk* g2;
      e3 += exk* g3;
      e4 += exk* g4;
      e5 += exs* g1;
      e6 += exs* g2;
      e7 += exs* g3;
      e8 += exs* g4;
      e9 += exk* f1+ exs*f2;
    }
  } /* for( i1 = 0; i1 < nint; i1++ ) */
  
  e[0]= e1;
  e[1]= e2;
  e[2]= e3;
  e[3]= e4;
  e[4]= e5;
  e[5]= e6;
  e[6]= e7;
  e[7]= e8;
  e[8]= e9;
  xj= xxj;
  yj= xyj;
  zj= xzj;
  m_s= xs;
}

/*-----------------------------------------------------------------------*/

/* impedance_print sets up the print formats for impedance loading */
void nec_context::impedance_print( int in1, int in2, int in3, nec_float fl1, nec_float fl2,
    nec_float fl3, nec_float fl4, nec_float fl5, nec_float fl6, const char *ia)
{
  /* record to be output and buffer used to make it */
  std::string record;
  char buf[15];
  int in[3], i1, i;
  nec_float fl[6];
  
  in[0]= in1;
  in[1]= in2;
  in[2]= in3;
  fl[0]= fl1;
  fl[1]= fl2;
  fl[2]= fl3;
  fl[3]= fl4;
  fl[4]= fl5;
  fl[5]= fl6;
  
  /* integer format */
  i1=0;
  record = "\n ";
  
  if ( (in1 == 0) && (in2 == 0) && (in3 == 0) )
  {
    record += " ALL";
    i1=1;
  }
  
  for( i = i1; i < 3; i++ )
  {
    if ( in[i] == 0)
      record += "     ";
    else
    {
      sprintf( buf, "%5d", in[i] );
      record += buf;
    }
  }
  
  /* floating point format */
  for( i = 0; i < 6; i++ )
  {
    if ( fabs( fl[i]) >= 1.0e-20 )
    {
      sprintf( buf, " %11.4E", fl[i] );
      record += buf;
    }
    else
      record += "            ";
  }
  
  record += "   ";
  record += ia;
  m_output.string(record.c_str() );
}

/*-----------------------------------------------------------------------*/

/* fill incident field array for charge discontinuity voltage source */
void nec_context::qdsrc( int is, nec_complex v, complex_array& e )
{
  static nec_complex s_CCJ(0.0,-0.01666666667);
  
  int i, jx, j, jp1, ipr, i1;
  nec_float xi, yi, zi, ai, cabi, sabi, salpi, tx, ty, tz;
  nec_complex curd, etk, ets, etc;
  
  is--;
  i= m_geometry->icon1[is];
  m_geometry->icon1[is]=0;
  m_geometry->tbf( is+1,0);
  m_geometry->icon1[is]= i;
  m_s = m_geometry->segment_length[is]*.5;
  curd= s_CCJ * v/(( log(2.* m_s/ m_geometry->segment_radius[is])-1.)*( m_geometry->bx[m_geometry->jsno-1]*
    cos( two_pi() * m_s)+ m_geometry->cx[m_geometry->jsno-1]* sin( two_pi() * m_s))* _wavelength);
  vqds[nqds]= v;
  iqds[nqds]= is+1;
  nqds++;
  
  for( jx = 0; jx < m_geometry->jsno; jx++ ) {
    j= m_geometry->jco[jx]-1;
    jp1 = j+1;
    m_s= m_geometry->segment_length[j];
    m_b= m_geometry->segment_radius[j];
    xj= m_geometry->x[j];
    yj= m_geometry->y[j];
    zj= m_geometry->z[j];
    cabj= m_geometry->cab[j];
    sabj= m_geometry->sab[j];
    salpj= m_geometry->salp[j];
  
    if ( m_use_exk == true) {
    ipr= m_geometry->icon1[j];
  
    if (ipr > PCHCON)
      ind1=2;
    else if ( ipr < 0 ) {
      ipr = - ipr;
      ipr--;
      if ( -m_geometry->icon1[ipr-1] != jp1 ) {
        ind1=2;
      } else {
        ind1 = m_geometry->test_ek_approximation(j,ipr);
      }
    }  /* if ( ipr < 0 ) */
    else if ( ipr == 0 )
      ind1=1;
    else /* ipr > 0 */
    {
      ipr--;
      if ( ipr != j ) {
        if ( m_geometry->icon2[ipr] != jp1) {
          ind1=2;
        } else {
          ind1 = m_geometry->test_ek_approximation(j,ipr);
        }
      } /* if ( ipr != j ) */
      else {
        if ( (cabj* cabj+ sabj* sabj) > 1.0e-8)
          ind1=2;
        else
          ind1=0;
      }
    } /* else */
  
    ipr= m_geometry->icon2[j];
    if (ipr > PCHCON) ind2=2;
    else if ( ipr < 0 ) {
      ipr = -ipr;
      ipr--;
      if ( -m_geometry->icon2[ipr] != jp1 ) {
        ind1=2;
      } else {
        ind1 = m_geometry->test_ek_approximation(j,ipr);
      }
    } /* if ( ipr < 0 ) */
    else if ( ipr == 0 )
      ind2=1;
    else /* ipr > 0 */
    {
      ipr--;
      if ( ipr != j ) {
        if ( m_geometry->icon1[ipr] != jp1) {
          ind2=2;
        } else {
          ind2 = m_geometry->test_ek_approximation(j,ipr);
        }
      } /* if ( ipr != j )*/
      else {
        if ( (cabj* cabj+ sabj* sabj) > 1.0e-8)
          ind1=2;
        else
          ind1=0;
      }
    } /* else */
  
    } /* if ( m_use_exk == true) */
  
    int n = m_geometry->n_segments;
    for( i = 0; i < n; i++ ) {
      xi= m_geometry->x[i];
      yi= m_geometry->y[i];
      zi= m_geometry->z[i];
      ai= m_geometry->segment_radius[i];
      efld( xi, yi, zi, ai, (i != j));
      cabi= m_geometry->cab[i];
      sabi= m_geometry->sab[i];
      salpi= m_geometry->salp[i];
      etk= exk* cabi+ eyk* sabi+ ezk* salpi;
      ets= exs* cabi+ eys* sabi+ ezs* salpi;
      etc= exc* cabi+ eyc* sabi+ ezc* salpi;
      e[i]= e[i]-( etk* m_geometry->ax[jx]+ ets* m_geometry->bx[jx]+ etc* m_geometry->cx[jx])* curd;
    }
  
    int m = m_geometry->m;
    if ( m != 0) {
      i1= n-1;
      for( i = 0; i < m; i++ ) {
        xi= m_geometry->px[i];
        yi= m_geometry->py[i];
        zi= m_geometry->pz[i];
        hsfld( xi, yi, zi,0.);
        i1++;
        tx= m_geometry->t2x[i];
        ty= m_geometry->t2y[i];
        tz= m_geometry->t2z[i];
        etk= exk* tx+ eyk* ty+ ezk* tz;
        ets= exs* tx+ eys* ty+ ezs* tz;
        etc= exc* tx+ eyc* ty+ ezc* tz;
        e[i1] += ( etk* m_geometry->ax[jx]+ ets* m_geometry->bx[jx]+ etc* m_geometry->cx[jx] )* curd* m_geometry->psalp[i];
        i1++;
        tx= m_geometry->t1x[i];
        ty= m_geometry->t1y[i];
        tz= m_geometry->t1z[i];
        etk= exk* tx+ eyk* ty+ ezk* tz;
        ets= exs* tx+ eys* ty+ ezs* tz;
        etc= exc* tx+ eyc* ty+ ezc* tz;
        e[i1] += ( etk* m_geometry->ax[jx]+ ets* m_geometry->bx[jx]+ etc* m_geometry->cx[jx])* curd* m_geometry->psalp[i];
      }  
    } /* if ( m != 0) */
  
    if ( nload > 0 )
      e[j] += zarray[j]* curd*(m_geometry->ax[jx]+ m_geometry->cx[jx]);
  
  } /* for( jx = 0; jx < m_geometry->jsno; jx++ ) */
}
  
/*-----------------------------------------------------------------------*/

#if 1

/*
Results: (epsilon = 81)
cond  --------- NEC-4D ------------- ------------- NEC-2D -----------
S/m  Rin     Xin    Avg. Gain  Rin     Xin     Avg. Gain

0.001  7.69123E+01  1.26695E+02  0.160    7.90906E+01  1.29476E+02  0.155
0.003  7.21069E+01  1.25749E+02  0.167    7.38409E+01   1.28334E+02   0.162
0.01  6.06862E+01  1.19449E+02  0.218    6.02083E+01   1.22146E+02  0.202
0.03  4.92272E+01  1.02255E+02  0.289    4.36457E+01  1.07793E+02   0.323
0.1  4.10716E+01  7.86277E+01  0.409    2.84842E+01   9.36058E+01   0.583
0.3  2.78348E+01  6.67538E+01  0.664    1.96676E+01   8.87374E+01   0.929
1.0  -4.87995E+00  1.01535E+02  -4.060    1.59479E+01   8.70547E+01   1.218
5.0  -4.21415E+01  3.45929E+02  -0.492    .42931E+01   8.59620E+01   1.472
*/
void nec_context::rom2( nec_float a, nec_float b, complex_array& sum, nec_float dmin )
{
/*
      Z=A
      ZE=B
      S=B-A
      IF (S.GE.0.) GO TO 1
      WRITE(3,18)
      STOP
*/
  int nts = 4, nx = 1;
  const int n = 9;
  nec_float rx = 1e-4;
  nec_float dz = 0.0;
  nec_float dzot= dz*.5;

  nec_float _z = a;
  nec_float ze = b;
  nec_float _s = b - a;
  if ( _s < 0.0)
  {
    throw new nec_exception("ERROR - B LESS THAN A IN ROM2");
  }


/*
1     EP=S/(1.E4*NM)
      ZEND=ZE-EP
      DO 2 I=1,N
2     SUM(I)=(0.,0.)
      NS=NX
      NT=0
      CALL SFLDS (Z,G1)
*/
  complex_array g1(9), g2(9), g3(9), g4(9), g5(9);
  complex_array t01(9), t10(9), t20(9), t11(9);

  /*! tolerance for hitting upper limit */
  const nec_float ep = _s/(1.0e4 * m_geometry->n_plus_m());

  /*! upper limit */
  const nec_float zend = ze - ep;
  
  vector_fill(sum,0,n,cplx_00());
  
  int ns = nx;
  /*! counter to control increasing subinterval size */
  int nt = 0;
  sflds( _z, g1);

/*
3     DZ=S/NS
      IF (Z+DZ.LE.ZE) GO TO 4
      DZ=ZE-Z
      IF (DZ.LE.EP) GO TO 17
4     DZOT=DZ*.5
      CALL SFLDS (Z+DZOT,G3)
      CALL SFLDS (Z+DZ,G5)
*/
loop3:  
  
  dz = _s / ns;
  if ( _z + dz >= ze)
  {
    dz = ze - _z;
    if ( dz <= ep)
      return;
  }

  dzot= dz*.5;
  sflds( _z + dzot, g3);
  sflds( _z + dz, g5);
  
  
/*
5     TMAG1=0.
      TMAG2=0.
C
C     EVALUATE 3 POINT ROMBERG RESULT AND TEST CONVERGENCE.
C
      DO 6 I=1,N
      T00=(G1(I)+G5(I))*DZOT
      T01(I)=(T00+DZ*G3(I))*.5
      T10(I)=(4.*T01(I)-T00)/3.
      IF (I.GT.3) GO TO 6
      TR=DREAL(T01(I))
      TI=DIMAG(T01(I))
      TMAG1=TMAG1+TR*TR+TI*TI
      TR=DREAL(T10(I))
      TI=DIMAG(T10(I))
      TMAG2=TMAG2+TR*TR+TI*TI
6     CONTINUE
*/
loop5:
  for (int i=0;i<n;i++)
  {
    nec_complex t00 = (g1[i]+ g5[i]) * dzot;
    t01[i] = (t00 + dz * g3[i]) * 0.5;
    t10[i] = (4.0 * t01[i] - t00)/3.;
  }

/*
      TMAG1=SQRT(TMAG1)
      TMAG2=SQRT(TMAG2)
      CALL TEST(TMAG1,TMAG2,TR,0.D0,0.D0,TI,DMIN)
*/
  nec_float tmag1 = sqrt( norm(t01[0]) + norm(t01[1]) + norm(t01[2]) );
  nec_float tmag2 = sqrt( norm(t10[0]) + norm(t10[1]) + norm(t10[2]) );

  nec_float tr = test_simple( tmag1, tmag2, dmin);

/*
      IF(TR.GT.RX)GO TO 8
      DO 7 I=1,N
7     SUM(I)=SUM(I)+T10(I)
      NT=NT+2
      GO TO 12
*/
  if (tr <= rx)
  {
    for (int i = 0; i < n; i++ )
      sum[i] += t10[i];
    nt += 2;
    goto loop12;
  }

/*
8     CALL SFLDS (Z+DZ*.25,G2)
      CALL SFLDS (Z+DZ*.75,G4)
      TMAG1=0.
      TMAG2=0.
C
C     EVALUATE 5 POINT ROMBERG RESULT AND TEST CONVERGENCE.
C
      DO 9 I=1,N
      T02=(T01(I)+DZOT*(G2(I)+G4(I)))*.5
      T11=(4.*T02-T01(I))/3.
      T20(I)=(16.*T11-T10(I))/15.
      IF (I.GT.3) GO TO 9
      TR=DREAL(T11)
      TI=DIMAG(T11)
      TMAG1=TMAG1+TR*TR+TI*TI
      TR=DREAL(T20(I))
      TI=DIMAG(T20(I))
      TMAG2=TMAG2+TR*TR+TI*TI
9     CONTINUE
      TMAG1=SQRT(TMAG1)
      TMAG2=SQRT(TMAG2)
      CALL TEST(TMAG1,TMAG2,TR,0.D0,0.D0,TI,DMIN)
*/
  sflds( _z+ dz*.25, g2);
  sflds( _z+ dz*.75, g4);
  
  /* Evaluate 5 point Romberg result and test convergence. */
  for (int i = 0; i < n; i++ )
  {
    nec_complex t02 = (t01[i]+ dzot*( g2[i]+ g4[i]))*0.5;
    t11[i] = (4.0 * t02- t01[i] )/3.0;
    t20[i] = (16.* t11[i] - t10[i])/15.0;
  }

  tmag1 = sqrt( norm(t11[0]) + norm(t11[1]) + norm(t11[2]) );
  tmag2 = sqrt( norm(t20[0]) + norm(t20[1]) + norm(t20[2]) );

  tr = test_simple( tmag1, tmag2, dmin);
/*
      IF(TR.GT.RX)GO TO 14
10    DO 11 I=1,N
11    SUM(I)=SUM(I)+T20(I)
      NT=NT+1
12    Z=Z+DZ
      IF (Z.GT.ZEND) GO TO 17
      DO 13 I=1,N
13    G1(I)=G5(I)
      IF (NT.LT.NTS.OR.NS.LE.NX) GO TO 3
      NS=NS/2
      NT=1
      GO TO 3
*/
  if (tr <= rx)
  {
loop10:
    for (int i = 0; i < n; i++ )
      sum[i] += t20[i];

    nt = nt+1;
  
loop12:
    _z = _z + dz;
    if ( _z > zend)
      return;
  
    for (int i = 0; i < n; i++ )
      g1[i] = g5[i];

    if ( (nt >= nts) && (ns > nx) )
    {
      /* Double step size */
      ns = ns/2;
      nt = 1;
    }
    goto loop3;
  }

/*
14    NT=0
      IF (NS.LT.NM) GO TO 15
      WRITE(3,19)  Z
      GO TO 10
15    NS=NS*2
      DZ=S/NS
      DZOT=DZ*.5
      DO 16 I=1,N
      G5(I)=G3(I)
16    G3(I)=G2(I)
      GO TO 5
*/
  nt=0;
  if ( ns >= m_geometry->n_plus_m() )
  {
    nec_error_mode em(m_output);
    m_output.string("ROM2 -- STEP SIZE LIMITED AT Z = ");
    m_output.real_out(12,5,_z);
    m_output.endl();
    goto loop10;
  }

  ns = ns*2;
  dz = _s/ ns;
  dzot = dz*.5;
  
  for (int i = 0; i < n; i++ )
  {
    g5[i] = g3[i];
    g3[i] = g2[i];
  }
  goto loop5;
}
/*
17    CONTINUE
      RETURN
C
18    FORMAT (30H ERROR - B LESS THAN A IN ROM2)
19    FORMAT (33H ROM2 -- STEP SIZE LIMITED AT Z =,1P,E12.5)
      END
*/

#else

/*!
  For the Sommerfeld ground option, rom2 integrates over the source 
  segment to obtain the total field due to ground.  The method of
  variable interval width Romberg integration is used.  There are 9
  field components - the x, y, and z components due to constant, 
  sine, and cosine current distributions.
  
  \param a lower limit of integral
  \param b upper limit of integral
  \param dmin Set in EFLD to 1% of the magnitude of the electric field
*/
void nec_context::rom2( nec_float a, nec_float b, complex_array& sum, nec_float dmin )
{
  ASSERT(sum.size() == 9);
  
  bool recalculate_fields = true;
  static bool step_warning_issued = false;
  
  int nts = 4, nx = 1, n = 9;

  /*! subinterval size */
  nec_float dz=0.0;

  /*! dz divided by 2.0 */
  nec_float dzot=0.0;

  /*! relative error limit */
  nec_float rx = 1.0e-8;
  
  complex_array g1(9), g2(9), g3(9), g4(9), g5(9);
  complex_array t01(9), t10(9), t20(9);
  
  nec_float _z = a;
  const nec_float ze = b;
  const nec_float _s = b - a;
  
  if ( _s < 0.0)
  {
    throw new nec_exception("ERROR - B LESS THAN A IN ROM2");
  }
  
  /*! tolerance for hitting upper limit */
  const nec_float ep = _s/(1.0e4 * m_geometry->n_plus_m());

  /*! upper limit */
  const nec_float zend = ze - ep;
  
  vector_fill(sum,0,n,cplx_00());
  
  int ns = nx;

  /*! counter to control increasing subinterval size */
  int nt = 0;
  sflds( _z, g1);

  while ( true )
  {
    if ( recalculate_fields )
    {
      // Now we do the same housekeeping before looping back
      _z = _z + dz;
      if ( _z > zend)
        return;
    
      for (int i = 0; i < n; i++ )
        g1[i] = g5[i];
  
      if ( (nt >= nts) && (ns > nx) )
      {
        /* Double step size */
        ns = ns/2;
        nt = 1;
      }

      dz = _s / ns;
      if ( _z + dz >= ze)
      {
        dz = ze - _z;
        if ( dz <= ep)
          return;
      }
    
      dzot= dz*.5;
      sflds( _z + dzot, g3);
      sflds( _z + dz, g5);
    }
  
    // Evaluate 3-point Romberg result and test convergence. 
    for (int i = 0; i < n; i++ )
    {
      nec_complex t00 = (g1[i]+ g5[i]) * dzot;
      t01[i] = (t00 + dz * g3[i]) * 0.5;
      t10[i] = (4.0 * t01[i] - t00)/3.;
    }
    
    nec_float tmag1 = sqrt( norm(t01[0]) + norm(t01[1]) + norm(t01[2]) );
    nec_float tmag2 = sqrt( norm(t10[0]) + norm(t10[1]) + norm(t10[2]) );

    nec_float tr = test_simple( tmag1, tmag2, dmin);
    
    if ( tr <= rx)
    {
      nt += 2;

      for (int i = 0; i < n; i++ )
        sum[i] += t10[i];
    
      recalculate_fields = true;
      continue;
    
    } /* if ( tr <= rx) */
  
    sflds( _z+ dz*.25, g2);
    sflds( _z+ dz*.75, g4);
    
    tmag1 = 0.0;
    tmag2 = 0.0;
  
    /* Evaluate 5 point Romberg result and test convergence. */
    for (int i = 0; i < n; i++ )
    {
      nec_complex t02 = (t01[i]+ dzot*( g2[i]+ g4[i]))*0.5;
      nec_complex t11 = (4.0 * t02- t01[i] )/3.0;
      t20[i] = (16.* t11- t10[i])/15.0;
      
      if (i <= 2)
      {
        tmag1 += norm(t11);
        tmag2 += norm(t20[i]);
      }
    }
  
    tmag1 = sqrt(tmag1);
    tmag2 = sqrt(tmag2);
    tr = test_simple( tmag1, tmag2, dmin);

    if ( tr > rx)
    {
      nt=0;
      if ( ns <= m_geometry->n_plus_m() )
      {
        // halve step size
        ns = ns*2;
        dz = _s/ ns;
        dzot = dz*.5;
      
        for (int i = 0; i < n; i++ )
        {
          g5[i] = g3[i];
          g3[i] = g2[i];
        }
      
        recalculate_fields=false;
        continue;
      }
    
      nec_error_mode em(m_output);
      m_output.string("ROM2 -- STEP SIZE LIMITED AT Z = ");
      m_output.real_out(12,5,_z);
      m_output.endl();
      
      if (false == step_warning_issued)
      {
        m_output.line("About the above warning:");
        m_output.line("Probably caused by a wire too close to the ground in the Somerfeld/");
        m_output.line("Norton ground method.  Execution continues but results may be inaccurate.");
        step_warning_issued = true;
      }
    } /* if ( tr <= rx) */
  
    nt = nt+1;

    for (int i = 0; i < n; i++ )
      sum[i] += t20[i];
  
    recalculate_fields = true;
    
  } /* while( true ) */
}
#endif

/*-----------------------------------------------------------------------*/

/*
  sflds returns the field due to ground for a current element on
  the source segment at t relative to the segment center.
*/
void nec_context::sflds(const nec_float t, complex_array& e )
{
  static nec_complex __const1(0.0,4.771341189);

  nec_float xt, yt, zt, rhx, rhy, rhs, rho, phx, phy;
  nec_float cph, sph, zphs, r2s, rk, sfac, thet;
  nec_complex  erv, ezv, erh, ezh, eph;

  xt= xj + t* cabj;
  yt= yj + t* sabj;
  zt= zj + t* salpj;
  rhx= xo- xt;
  rhy= yo- yt;
  rhs= rhx* rhx+ rhy* rhy;
  rho= sqrt( rhs);

  if ( rho <= 0.0) {
    rhx=1.0;
    rhy=0.0;
    phx=0.0;
    phy=1.0;
  } else {
    rhx= rhx/ rho;
    rhy= rhy/ rho;
    phx= -rhy;
    phy= rhx;
  }

  cph= rhx* xsn+ rhy* ysn;
  sph= rhy* xsn- rhx* ysn;

  if ( fabs( cph) < 1.0e-10)
    cph=0.0;
  if ( fabs( sph) < 1.0e-10)
    sph=0.0;

  ground_wave.zph = zo+ zt;
  zphs= ground_wave.zph* ground_wave.zph;
  r2s= rhs+ zphs;
  ground_wave.r2= sqrt( r2s);
  rk= ground_wave.r2* two_pi();
  ground_wave.xx2 = nec_complex( cos( rk),-sin( rk));

  /*  Use Norton approximation for field due to ground.  Current is
    lumped at segment center with current moment for constant, sine,
    or cosine distribution. */
  if ( isnor != 1) {
    ground_wave.zmh=1.0;
    ground_wave.r1=1.;
    ground_wave.xx1=0.;
    gwave(erv, ezv, erh, ezh, eph, ground_wave);

    nec_complex et = -__const1 * ground.frati* ground_wave.xx2/( r2s* ground_wave.r2);
    nec_complex er = 2.* et* nec_complex(1.0, rk);
    et= et* nec_complex(1.0 - rk* rk, rk);
    nec_complex hrv = ( er+ et)* rho* ground_wave.zph/ r2s;
    nec_complex hzv = ( zphs* er- rhs* et)/ r2s;
    nec_complex hrh = ( rhs* er- zphs* et)/ r2s;
    erv= erv- hrv;
    ezv= ezv- hzv;
    erh= erh+ hrh;
    ezh= ezh+ hrv;
    eph= eph+ et;
    erv= erv* salpj;
    ezv= ezv* salpj;
    erh= erh* sn* cph;
    ezh= ezh* sn* cph;
    eph= eph* sn* sph;
    erh= erv+ erh;
    e[0]=( erh* rhx+ eph* phx)* m_s;
    e[1]=( erh* rhy+ eph* phy)* m_s;
    e[2]=( ezv+ ezh)* m_s;
    e[3]=0.;
    e[4]=0.;
    e[5]=0.;
    sfac= pi()* m_s;
    sfac= sin( sfac)/ sfac;
    e[6]= e[0]* sfac;
    e[7]= e[1]* sfac;
    e[8]= e[2]* sfac;

    return;
  } /* if ( isnor != 1) */

  /* Interpolate in Sommerfeld field tables */
  if ( rho >= 1.0e-12)
    thet= atan( ground_wave.zph/ rho);
  else
    thet= pi_two();

  /*  Combine vertical and horizontal components and convert
    to x,y,z components. multiply by exp(-jkr)/r.
  */
  ground.ggrid_interpolate( ground_wave.r2, thet, &erv, &ezv, &erh, &eph );
  ground_wave.xx2= ground_wave.xx2 / ground_wave.r2;
  sfac= sn* cph;
  erh= ground_wave.xx2*( salpj* erv+ sfac* erh);
  ezh= ground_wave.xx2*( salpj* ezv- sfac* erv);
  /* x,y,z fields for constant current */
  eph= sn* sph* ground_wave.xx2* eph;
  e[0]= erh* rhx+ eph* phx;
  e[1]= erh* rhy+ eph* phy;
  e[2]= ezh;
  /* x,y,z fields for sine current */
  rk= two_pi() * t;
  sfac= sin( rk);
  e[3]= e[0]* sfac;
  e[4]= e[1]* sfac;
  e[5]= e[2]* sfac;
  /* x,y,z fields for cosine current */
  sfac= cos( rk);
  e[6]= e[0]* sfac;
  e[7]= e[1]* sfac;
  e[8]= e[2]* sfac;
}
/*-----------------------------------------------------------------------*/

/*!\brief Calculates the electric field due to unit current in the t1 and t2 directions on a patch
\param ground_reflection If true, then calculate the field reflected from the ground. (was ipgnd == 2)
*/
void nec_context::unere( nec_float xob, nec_float yob, nec_float zob, bool ground_reflection )
{
  nec_float px, py, cth;
  nec_complex rrv, rrh, edp;
  
  nec_float zr = zj;
  nec_float t1zr = t1zj;
  nec_float t2zr = t2zj;
  
  if ( ground_reflection)
  {
    zr = - zr;
    t1zr = -t1zr;
    t2zr = -t2zr;
  }
  
  nec_float rx = xob- xj;
  nec_float ry = yob- yj;
  nec_float rz = zob- zr;
  nec_float r2 = rx*rx + ry*ry + rz*rz;
  
  if ( r2 <= 1.0e-20)
  {
    exk=cplx_00();
    eyk=cplx_00();
    ezk=cplx_00();
    exs=cplx_00();
    eys=cplx_00();
    ezs=cplx_00();
    return;
  }
  
  nec_float r = sqrt(r2);
  nec_float tt1 = -two_pi() * r;
  nec_float tt2 = tt1 * tt1;
  nec_float rt = r2*r;
  nec_complex er = nec_complex( sin(tt1),-cos(tt1))*( CONST2* m_s);
  nec_complex q1= nec_complex(tt2 - 1.0, tt1)* er/ rt;
  nec_complex q2= nec_complex(3.0- tt2,-3.0*tt1)* er/( rt* r2);
  er = q2*( t1xj* rx+ t1yj* ry+ t1zr* rz);
  exk= q1* t1xj+ er* rx;
  eyk= q1* t1yj+ er* ry;
  ezk= q1* t1zr+ er* rz;
  er= q2*( t2xj* rx+ t2yj* ry+ t2zr* rz);
  exs= q1* t2xj+ er* rx;
  eys= q1* t2yj+ er* ry;
  ezs= q1* t2zr+ er* rz;
  
  if ( !ground_reflection)
    return;
  
  // handle the ground_reflection
  if (  ground.type_perfect() ) // (ground.iperf == 1)
  {
    exk = - exk;
    eyk = - eyk;  
    
    ezk = - ezk;
    exs = - exs;
    eys = - eys;
    ezs = - ezs;
    return;
  }
  
  nec_float xymag = sqrt(rx*rx + ry*ry);
  if ( xymag <= 1.0e-6)
  {
    px=0.;
    py=0.;
    cth=1.;
    rrv=cplx_10();
  }
  else
  {
    px = - ry/ xymag;
    py= rx/ xymag;
    cth= rz/ sqrt( xymag* xymag+ rz* rz);
    rrv= sqrt(1.0- ground.get_zrati_sqr() * (1.0 - cth*cth));
  }
  
  rrh= ground.zrati* cth;
  rrh=( rrh- rrv)/( rrh+ rrv);
  rrv= ground.zrati* rrv;
  rrv = -( cth- rrv)/( cth+ rrv);
  
  edp=( exk* px+ eyk* py)*( rrh- rrv);
  exk= exk* rrv+ edp* px;
  eyk= eyk* rrv+ edp* py;
  ezk= ezk* rrv;
  edp=( exs* px+ eys* py)*( rrh- rrv);
  exs= exs* rrv+ edp* px;
  eys= eys* rrv+ edp* py;
  ezs= ezs* rrv;
}

/*-----------------------------------------------------------------------*/


/* zint computes the internal impedance of a circular wire */
nec_complex zint_old( nec_float sigl, nec_float rolam );
nec_complex zint_old( nec_float sigl, nec_float rolam )
{
#define cc1  nec_complex( 6.0e-7,     + 1.9e-6)
#define cc2  nec_complex(-3.4e-6,     + 5.1e-6)
#define cc3  nec_complex(-2.52e-5,    + 0.0)
#define cc4  nec_complex(-9.06e-5 ,   - 9.01e-5)
#define cc5  nec_complex( 0.,         - 9.765e-4)
#define cc6  nec_complex(.0110486,    - .0110485)
#define cc7  nec_complex( 0.,         - .3926991)
#define cc8  nec_complex( 1.6e-6,     - 3.2e-6)
#define cc9  nec_complex( 1.17e-5,    - 2.4e-6)
#define cc10  nec_complex( 3.46e-5,    + 3.38e-5)
#define cc11  nec_complex( 5.0e-7,     + 2.452e-4)
#define cc12  nec_complex(-1.3813e-3,  + 1.3811e-3)
#define cc13  nec_complex(-6.25001e-2, - 1.0e-7)
#define cc14  nec_complex(.7071068,    + .7071068)
#define cn  cc14

#define th(d) ( (((((cc1*(d)+cc2)*(d)+cc3)*(d)+cc4)*(d)+cc5)*(d)+cc6)*(d) + cc7 )
#define ph(d) ( (((((cc8*(d)+cc9)*(d)+cc10)*(d)+cc11)*(d)+cc12)*(d)+cc13)*(d)+cc14 )
#define f(d)  ( sqrt(pi_two()/(d))*exp(-cn*(d)+th(-8./x)) )
#define g(d)  ( exp(cn*(d)+th(8./x))/sqrt(two_pi()*(d)) )

  nec_complex br1, br2, zint;
  nec_float x, y, s, ber, bei;
  nec_float tpcmu = 2.368705e+3;
  nec_float cmotp = 60.00;

  x= sqrt( tpcmu* sigl)* rolam;
  if ( x <= 110.)
  {
    if ( x <= 8.)
    {
      y= x/8.;
      y= y* y;
      s= y* y;

      ber=((((((-9.01e-6* s+1.22552e-3)* s-.08349609)* s+ 2.6419140)*
        s-32.363456)* s+113.77778)* s-64.)* s+1.;

      bei=((((((1.1346e-4* s-.01103667)* s+.52185615)* s-10.567658)*
        s+72.817777)* s-113.77778)* s+16.)* y;

      br1= nec_complex( ber, bei);

      ber=(((((((-3.94e-6* s+4.5957e-4)* s-.02609253)* s+ .66047849)*
    s-6.0681481)* s+14.222222)* s-4.)* y)* x;

      bei=((((((4.609e-5* s-3.79386e-3)* s+.14677204)* s- 2.3116751)*
        s+11.377778)* s-10.666667)* s+.5)* x;

      br2= nec_complex( ber, bei);
      br1= br1/ br2;
      zint= cplx_01()* sqrt( cmotp/sigl )* br1/ rolam;

      return( zint );

    } // if ( x <= 8.) 

    br2= cplx_01()* f(x)/ pi();
    br1= g( x)+ br2;
    br2= g( x)* ph(8./ x)- br2* ph(-8./ x);
    br1= br1/ br2;
    zint= cplx_01()* sqrt( cmotp/ sigl)* br1/ rolam;

    return( zint );

  } // if ( x <= 110.) 

  br1= nec_complex(.70710678,-.70710678);
  zint= cplx_01()* sqrt( cmotp/ sigl)* br1/ rolam;

  return( zint );
}


/* zint computes the internal impedance of a circular wire */
nec_complex nec_context::zint( nec_float sigl, nec_float rolam )
{
#define cc1  nec_complex( 6.0e-7,     + 1.9e-6)
#define cc2  nec_complex(-3.4e-6,     + 5.1e-6)
#define cc3  nec_complex(-2.52e-5,    + 0.0)
#define cc4  nec_complex(-9.06e-5 ,   - 9.01e-5)
#define cc5  nec_complex( 0.,         - 9.765e-4)
#define cc6  nec_complex(.0110486,    - .0110485)
#define cc7  nec_complex( 0.,         - .3926991)
#define cc8  nec_complex( 1.6e-6,     - 3.2e-6)
#define cc9  nec_complex( 1.17e-5,    - 2.4e-6)
#define cc10  nec_complex( 3.46e-5,    + 3.38e-5)
#define cc11  nec_complex( 5.0e-7,     + 2.452e-4)
#define cc12  nec_complex(-1.3813e-3,  + 1.3811e-3)
#define cc13  nec_complex(-6.25001e-2, - 1.0e-7)
#define cc14  nec_complex(.7071068,    + .7071068)
#define cn  cc14

#define th(d) ( (((((cc1*(d)+cc2)*(d)+cc3)*(d)+cc4)*(d)+cc5)*(d)+cc6)*(d) + cc7 )
#define ph(d) ( (((((cc8*(d)+cc9)*(d)+cc10)*(d)+cc11)*(d)+cc12)*(d)+cc13)*(d)+cc14 )
#define f(d)  ( sqrt(pi_two()/(d))*exp(-cn*(d)+th(-8./x)) )
#define g(d)  ( exp(cn*(d)+th(8./x))/sqrt(two_pi()*(d)) )

  static nec_float tpcmu = 2.368705e+3;
  static nec_float cmotp = 60.00;

  nec_float x = sqrt(tpcmu * sigl) * rolam;
  
  // if x is zero, then we have a zero radius!
  ASSERT(x != 0.0);
    
  if (x > 110.0)
  {
    nec_complex br1 = nec_complex(0.70710678, -0.70710678);
    nec_complex ret= cplx_01()* sqrt( cmotp/ sigl)* br1/ rolam;
    
    ASSERT(ret == zint_old(sigl, rolam));
    return ret;
  }
  
  if (x > 8.0)
  {
    nec_complex br2 = cplx_01()* f(x)/ pi();
    nec_complex temp_gx = g(x);
    nec_complex br1 = temp_gx + br2;
    
    br2 = temp_gx * ph(8.0/x) - br2 * ph(-8.0/x);
    br1 = br1/ br2;
    
    nec_complex ret = cplx_01()* sqrt( cmotp/ sigl)* br1/ rolam;
  
    ASSERT(ret == zint_old(sigl, rolam));
    return ret;
  }
  
   nec_float x8 = x / 8.0;
  nec_float y = x8*x8;
  nec_float s = y*y;

  nec_float ber=((((((-9.01e-6* s+1.22552e-3)* s-.08349609)* s+ 2.6419140)*
    s-32.363456)* s+113.77778)* s-64.)* s+1.;

  nec_float bei=((((((1.1346e-4* s-.01103667)* s+.52185615)* s-10.567658)*
    s+72.817777)* s-113.77778)* s+16.)* y;

  nec_complex br1= nec_complex( ber, bei);

  ber=(((((((-3.94e-6*s + 4.5957e-4)*s - 0.02609253)*s + 0.66047849)*
  s - 6.0681481)*s + 14.222222)*s - 4.0)* y)* x;

  bei=((((((4.609e-5* s-3.79386e-3)* s+.14677204)* s- 2.3116751)*
    s+11.377778)* s-10.666667)* s+.5)* x;

  nec_complex br2= nec_complex( ber, bei);
  
  br1= br1/ br2;
  nec_complex ret = cplx_01()* sqrt( cmotp/sigl )* br1/ rolam;

  ASSERT(ret == zint_old(sigl, rolam));
  return ret;
}


/*

  fblock( np + 2 mp, n+2m, m_geometry->n_plus_2m * (m_geometry->np+2*m_geometry->mp),  m_geometry->m_ipsym)
*/
/* fblock sets parameters for out-of-core */
/* solution for the primary matrix (a) */
void nec_context::fblock( int nrow, int ncol, int imax, int ipsym ) {
  int ka, kk;
  
  if ( (nrow*ncol) <= imax)  {
    npblk= nrow;
    nlast= nrow;
    
    if ( nrow == ncol)  {
      icase=1;
      return;
    }
    else
      icase=2;
  } /* if ( nrow*ncol <= imax) */
  
  if ( (nop*nrow) != ncol)  {
    nec_stop("SYMMETRY ERROR - NROW: %d NCOL: %d", nrow, ncol );
  }
  
  /* set up symmetry_array matrix for rotational symmetry. */
  if ( ipsym <= 0)  {
    nec_float phaz = two_pi()/nop;
    
    for(int i = 1; i < nop; i++ )  {
      for(int j= i; j < nop; j++ )  {
        nec_float arg = phaz * (nec_float)i * (nec_float)j;
        symmetry_array[i+j*nop]= nec_complex( cos( arg), sin( arg));
        symmetry_array[j+i*nop]= symmetry_array[i+j*nop];
      }
    }
    return;
  } /* if ( ipsym <= 0) */
  
  /* set up symmetry_array matrix for plane symmetry */
  kk=1;
  symmetry_array[0]=cplx_10();
  
  int k_power = 2;
  for( ka = 1; k_power != nop; ka++ )
    k_power *= 2;
  
  for(int k = 0; k < ka; k++ )  {
    for(int i = 0; i < kk; i++ )  {
      for(int j = 0; j < kk; j++ )  {
        nec_complex deter = symmetry_array[i+j*nop];
        symmetry_array[i+(j+kk)*nop] = deter;
        symmetry_array[i+kk+(j+kk)*nop] = - deter;
        symmetry_array[i+kk+j*nop] = deter;
      }
    }
    kk *= 2;  
  } /* for( k = 0; k < ka; k++ ) */
}


/*! \brief gfld computes the radiated field including ground wave.

\param space_only Compute only the space wave (was ksymp == 1)
*/
void nec_context::gfld(nec_float rho, nec_float phi, nec_float rz,
    nec_complex *eth, nec_complex *epi,
    nec_complex *erd, bool space_only, nec_float in_wavelength )
{
  int i, k;
  nec_float b, r, thet, arg, phx, phy, rx, ry, dx, dy, dz, rix, riy, rhs, rhp;
  nec_float rhx, rhy, calp, cbet, sbet, cph, sph, el, rfl, riz, thx, thy, thz;
  nec_float rxyz, rnx, rny, rnz, omega, sill, top, bot, a, too, boo, c, rr, ri;
  nec_complex cix, ciy, ciz, exa, erv;
  nec_complex ezv, erh, eph, ezh, ex, ey;
  
  r= sqrt( rho*rho+ rz*rz );
  if ( (space_only) || (abs(ground.zrati) > .5) || (r > 1.e5) )  {
    /* computation of space wave only */
    if ( rz >= 1.0e-20)
      thet= atan( rho/ rz);
    else
      thet= pi()*.5;
  
    ffld(thet, phi, eth, epi, in_wavelength);
    arg= -two_pi() * r;
    exa= nec_complex( cos( arg), sin( arg))/ r;
    *eth= *eth* exa;
    *epi= *epi* exa;
    *erd=cplx_00();
    return;
  } /* if ( (space_only) && (abs(gound.zrati) > .5) && (r > 1.e5) ) */
  
  /* computation of space and ground waves. */
  ground_wave.set_u(ground.zrati);
  phx = - sin( phi);
  phy= cos( phi);
  rx= rho* phy;
  ry = - rho* phx;
  cix=cplx_00();
  ciy=cplx_00();
  ciz=cplx_00();
  
  /* summation of field from individual segments */
  for( i = 0; i < m_geometry->n_segments; i++ )  {
    dx= m_geometry->cab[i];
    dy= m_geometry->sab[i];
    dz= m_geometry->salp[i];
    rix= rx- m_geometry->x[i];
    riy= ry- m_geometry->y[i];
    rhs= rix* rix+ riy* riy;
    rhp= sqrt( rhs);
  
    if ( rhp >= 1.0e-6)  {
      rhx= rix/ rhp;
      rhy= riy/ rhp;
    } else {
      rhx=1.;
      rhy=0.;
    }
  
    calp=1.- dz* dz;
    if ( calp >= 1.0e-6)  {
      calp= sqrt( calp);
      cbet= dx/ calp;
      sbet= dy/ calp;
      cph= rhx* cbet+ rhy* sbet;
      sph= rhy* cbet- rhx* sbet;
    } else {
      cph= rhx;
      sph= rhy;
    }
  
    el= pi()* m_geometry->segment_length[i];
    rfl = -1.;
  
    /* Integration of (current)*(phase factor) over segment and image for 
       constant, sine, and cosine current distributions */
    for( k = 0; k < 2; k++ )  {
      rfl = - rfl;
      riz= rz- m_geometry->z[i]* rfl;
      rxyz= sqrt( rix* rix+ riy* riy+ riz* riz);
      rnx= rix/ rxyz;
      rny= riy/ rxyz;
      rnz= riz/ rxyz;
      omega = -( rnx* dx+ rny* dy+ rnz* dz* rfl);
      sill= omega* el;
      top= el+ sill;
      bot= el- sill;
    
      if ( fabs( omega) >= 1.0e-7)
        a=2.* sin( sill)/ omega;
      else
        a=(2.- omega* omega* el* el/3.)* el;
    
      if ( fabs( top) >= 1.0e-7)
        too= sin( top)/ top;
      else
        too=1.- top* top/6.;
    
      if ( fabs( bot) >= 1.0e-7)
        boo= sin( bot)/ bot;
      else
        boo=1.- bot* bot/6.;
    
      b= el*( boo- too);
      c= el*( boo+ too);
      rr= a* air[i]+ b* bii[i]+ c* cir[i];
      ri= a* aii[i]- b* bir[i]+ c* cii[i];
      arg= two_pi()*( m_geometry->x[i]* rnx+ m_geometry->y[i]* rny+ m_geometry->z[i]* rnz* rfl);
      exa= nec_complex( cos( arg), sin( arg))* nec_complex( rr, ri)/two_pi();
    
      if ( k != 1 )  {
        ground_wave.xx1= exa;
        ground_wave.r1= rxyz;
        ground_wave.zmh= riz;
        continue;
      }
    
      ground_wave.xx2 = exa;
      ground_wave.r2= rxyz;
      ground_wave.zph= riz;
    } /* for( k = 0; k < 2; k++ ) */
  
    /* call subroutine to compute the field */
    /* of segment including ground wave. */
    gwave(erv, ezv, erh, ezh, eph, ground_wave);
    erh= erh* cph* calp+ erv* dz;
    eph= eph* sph* calp;
    ezh= ezh* cph* calp+ ezv* dz;
    ex= erh* rhx- eph* rhy;
    ey= erh* rhy+ eph* rhx;
    cix= cix+ ex;
    ciy= ciy+ ey;
    ciz= ciz+ ezh;
  } /* for( i = 0; i < n; i++ ) */
  
  arg= -two_pi() * r;
  exa= nec_complex( cos( arg), sin( arg));
  cix= cix* exa;
  ciy= ciy* exa;
  ciz= ciz* exa;
  rnx= rx/ r;
  rny= ry/ r;
  rnz= rz/ r;
  thx= rnz* phy;
  thy = - rnz* phx;
  thz = - rho/ r;
  *eth= cix* thx+ ciy* thy+ ciz* thz;
  *epi= cix* phx+ ciy* phy;
  *erd= cix* rnx+ ciy* rny+ ciz* rnz;
}

/* ffld calculates the far zone radiated electric fields, */
/* the factor exp(j*k*r)/(r/lamda) not included */
void nec_context::ffld(nec_float thet, nec_float phi,
    nec_complex *eth, nec_complex *eph, nec_float in_wavelength )
{
  static nec_complex CONST3(0.0, -em::impedance() / four_pi()); // -29.97922085;

  int k, i;
  bool jump;
  nec_float phx, phy, roz, rozs, thx, thy, thz, rox, roy;
  nec_float tthet=0., darg=0., omega, el, sill, top, bot, a;
  nec_float too, boo, b, c, d, rr, ri, arg, dr;
  nec_complex cix, ciy, ciz, exa, ccx, ccy, ccz, cdp;
  nec_complex zrsin, rrv, rrh, rrv1, rrh1, rrv2, rrh2;
  nec_complex tix, tiy, tiz, ex, ey, ez;
  
  phx = - sin( phi);
  phy= cos( phi);
  roz= cos( thet);
  rozs= roz;
  thx= roz* phy;
  thy = - roz* phx;
  thz = - sin( thet);
  rox = - thz* phy;
  roy= thz* phx;
  
  jump = false;
  if ( m_geometry->n_segments != 0)  {
    /* loop for structure image if any */
    /* calculation of reflection coeffecients */
    for( k = 0; k < ground.ksymp; k++ )  {
    if ( k != 0 )  {
      /* for perfect ground */
      if (ground.type_perfect()) {
        rrv = -cplx_10();
        rrh = -cplx_10();
      } else {
        /* for infinite planar ground */
        zrsin= sqrt(1.- ground.get_zrati_sqr() * thz* thz);
        rrv = -( roz- ground.zrati * zrsin)/( roz+ ground.zrati* zrsin);
        rrh=( ground.zrati* roz- zrsin)/( ground.zrati* roz+ zrsin);
      }
    
      /* for the cliff problem, two reflction coefficients calculated */
      if ( ifar > 1)  {
        rrv1= rrv;
        rrh1= rrh;
        tthet= tan( thet);
      
        if ( ifar != 4)  {
          nec_complex zrati2 = ground.get_zrati2(in_wavelength);
          
          zrsin = sqrt(1.-  zrati2 *  zrati2 * thz* thz);
          rrv2 = -( roz-  zrati2* zrsin)/( roz+  zrati2* zrsin);
          rrh2 =(  zrati2* roz- zrsin)/(  zrati2* roz+ zrsin);
          darg = -two_pi() * 2.0 * ground.get_ch(in_wavelength) * roz;
        }
      } /* if ( ifar > 1) */
    
      roz = - roz;
      ccx= cix;
      ccy= ciy;
      ccz= ciz;
    } /* if ( k != 0 ) */
  
    cix=cplx_00();
    ciy=cplx_00();
    ciz=cplx_00();
  
    /* loop over structure segments */
    for( i = 0; i < m_geometry->n_segments; i++ )  {
      omega = -( rox* m_geometry->cab[i]+ roy* m_geometry->sab[i]+ roz* m_geometry->salp[i]);
      el= pi()* m_geometry->segment_length[i];
      sill= omega* el;
      top= el+ sill;
      bot= el- sill;
    
      if ( fabs( omega) >= 1.0e-7)
        a=2.* sin( sill)/ omega;
      else
        a=(2.- omega* omega* el* el/3.)* el;
    
      if ( fabs( top) >= 1.0e-7)
        too= sin( top)/ top;
      else
        too=1.- top* top/6.;
    
      if ( fabs( bot) >= 1.0e-7)
        boo= sin( bot)/ bot;
      else
        boo=1.- bot* bot/6.;
    
      b= el*( boo- too);
      c= el*( boo+ too);
      rr= a* air[i]+ b* bii[i]+ c* cir[i];
      ri= a* aii[i]- b* bir[i]+ c* cii[i];
      arg= two_pi()*( m_geometry->x[i]* rox+ m_geometry->y[i]* roy+ m_geometry->z[i]* roz);
    
      if ( (k != 1) || (ifar < 2) )  {
        /* summation for far field integral */
        exa= nec_complex( cos( arg), sin( arg))* nec_complex( rr, ri);
        cix= cix+ exa* m_geometry->cab[i];
        ciy= ciy+ exa* m_geometry->sab[i];
        ciz= ciz+ exa* m_geometry->salp[i];
        continue;
      }
    
      /* calculation of image contribution */
      /* in cliff and ground screen problems */
    
      /* specular point distance */
      dr= m_geometry->z[i]* tthet;
    
      d= dr* phy+ m_geometry->x[i];
      if ( ifar == 2)  {
        if (( ground.get_cl(in_wavelength) - d) > 0.0) {
          rrv= rrv1;
          rrh= rrh1;
        } else {
          rrv= rrv2;
          rrh= rrh2;
          arg= arg+ darg;
        }
      } else { /* if ( ifar == 2) */
        d= sqrt( d*d + (m_geometry->y[i]-dr*phx)*(m_geometry->y[i]-dr*phx) );
        if ( ifar == 3)  {
          if (( ground.get_cl(in_wavelength) - d) > 0.0)  {
            rrv= rrv1;
            rrh= rrh1;
          } else {
            rrv= rrv2;
            rrh= rrh2;
            arg= arg+ darg;
          }
        } else { /* if ( ifar == 3) */
          if (( ground.get_radial_wire_length_wavelengths() - d) >= 0.0)  {
            /* radial wire ground screen reflection coefficient */
            d += ground.t2;
            nec_complex zscrn= ground.m_t1 * d* log( d/ ground.t2);
            zscrn=( zscrn* ground.zrati)/( em::impedance() * ground.zrati+ zscrn);
            zrsin= sqrt(1.- zscrn* zscrn* thz* thz);
            rrv=( roz+ zscrn* zrsin)/(- roz+ zscrn* zrsin);
            rrh=( zscrn* roz+ zrsin)/( zscrn* roz- zrsin);
          } else  {  /* if (( ground.scrwl- d) < 0.) */
            if ( ifar == 4)  {
              rrv= rrv1;
              rrh= rrh1;
            } else  { /* if ( ifar == 4) */
              if ( ifar == 5)
                d= dr* phy+ m_geometry->x[i];
          
              if (( ground.get_cl(in_wavelength) - d) > 0.)  {
                rrv= rrv1;
                rrh= rrh1;
              } else {
                rrv= rrv2;
                rrh= rrh2;
                arg= arg+ darg;
              } /* if (( cl- d) > 0.) */
            } /* if ( ifar == 4) */
          } /* if (( ground.scrwl- d) < 0.) */
        } /* if ( ifar == 3) */
      } /* if ( ifar == 2) */
      
      /* contribution of each image segment modified by */
      /* reflection coef, for cliff and ground screen problems */
      exa= nec_complex( cos( arg), sin( arg))* nec_complex( rr, ri);
      tix= exa* m_geometry->cab[i];
      tiy= exa* m_geometry->sab[i];
      tiz= exa* m_geometry->salp[i];
      cdp=( tix* phx+ tiy* phy)*( rrh- rrv);
      cix= cix+ tix* rrv+ cdp* phx;
      ciy= ciy+ tiy* rrv+ cdp* phy;
      ciz= ciz- tiz* rrv;
    
      } /* for( i = 0; i < n; i++ ) */
  
    if ( k == 0 )
      continue;
  
    /* calculation of contribution of structure image for infinite ground */
    if ( ifar < 2) {
      cdp=( cix* phx+ ciy* phy)*( rrh- rrv);
      cix= ccx+ cix* rrv+ cdp* phx;
      ciy= ccy+ ciy* rrv+ cdp* phy;
      ciz= ccz- ciz* rrv;
    } else {
      cix= cix+ ccx;
      ciy= ciy+ ccy;
      ciz= ciz+ ccz;
    }
  
    } /* for( k=0; k < ground.ksymp; k++ ) */
  
    if ( m_geometry->m > 0) {
      jump = true;
    } else {
      *eth=( cix* thx+ ciy* thy+ ciz* thz)* CONST3;
      *eph=( cix* phx+ ciy* phy)* CONST3;
      return;
    }
  } /* if ( n != 0) */
  
  if ( ! jump )  {
    cix=cplx_00();
    ciy=cplx_00();
    ciz=cplx_00();
  }
  
  /* electric field components */
  roz= rozs;
  { // without ground 
    complex_array temp = current_vector.segment(m_geometry->n_segments, current_vector.size()-m_geometry->n_segments);
    m_geometry->fflds(rox, roy, roz, temp, &ex, &ey, &ez);
  }
  if (ground.present())  {
    // with ground
    nec_complex tempx, tempy, tempz;
    
    complex_array temp = current_vector.segment(m_geometry->n_segments, current_vector.size()-m_geometry->n_segments);
    m_geometry->fflds(rox, roy, -roz, temp, &tempx, &tempy, &tempz);
  
    if (ground.type_perfect())  {
      tempx = -tempx;
      tempy = -tempy;
      tempz = -tempz;
    } else {
      // get reflection co-efficient?? 
      rrv= sqrt(1.0 - ground.get_zrati_sqr() * thz* thz);
      rrh= ground.zrati * roz;
      rrh=( rrh- rrv)/( rrh+ rrv);
      rrv= ground.zrati * rrv;
      rrv = -( roz- rrv)/( roz+ rrv);
      
      *eth=( tempx* phx+ tempy* phy)*( rrh- rrv);
      tempx= tempx* rrv+ *eth* phx;
      tempy= tempy* rrv+ *eth* phy;
      tempz= tempz* rrv;
    } /* if (  ground.iperf == 1) */
  
    ex += tempx;
    ey += tempy;
    ez -= tempz;
  }
  
  ex += cix* CONST3;
  ey += ciy* CONST3;
  ez += ciz* CONST3;
  *eth= ex*thx + ey*thy + ez*thz;
  *eph= ex*phx + ey*phy;
}

