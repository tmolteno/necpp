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
#include "math_util.h"
#include "nec_exception.h"
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
#include "nec_ground.h"
#include "nec_output.h"
#include "misc.h" // for stop() function
#include "electromag.h"

/*
  Parse a GN card. The input parameters here are the fields of the
  GN card.
    
GN    NEAR GROUND, GROUND SCREEN, ADDED GROUND
  I1- -1=SET FREE SPACE (A), 0=REFL COEFF, 1=IDEAL (B), 2-SOMMERFIELD
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
void nec_ground::parse_gn(int itmp1, int itmp2,
  nec_float tmp1, nec_float tmp2,
  nec_float tmp3, nec_float tmp4,
  nec_float tmp5, nec_float tmp6
  )
{
  // iperf = -1 - nullifies ground parameters previously used and sets free- space condition.
  // The remainder of the card is left blank in this case. 
  
  if ( itmp1 == -1 ) { // nullify previous ground conditions.
    ksymp=1;
    radial_wire_count=0;
    iperf=0;
    return;
  }

  iperf = itmp1;
  ASSERT(iperf >= 0);
  
  radial_wire_count = itmp2;
  ksymp = 2;
  epsr = tmp1;
  sig = tmp2;

  if (radial_wire_count != 0) {
    if ( iperf == 2) {
      throw new nec_exception("RADIAL WIRE G.S. APPROXIMATION MAY NOT BE USED WITH SOMMERFELD GROUND OPTION");
    }

    radial_wire_length= tmp3;
    radial_wire_radius= tmp4;
    return;
  }

  if (0 == iperf) {
    if (epsr == 0.0) {
      throw new nec_exception("GROUND DIELECTRIC CONSTANT IS ZERO");
    }
    if (sig == 0.0) {
      throw new nec_exception("GROUND CONDUCTIVITY IS ZERO");
    }
  }
  
  setup_cliff(tmp3,tmp4,tmp5,tmp6);
}

/*
  Setup a cliff (two medium ground)
*/
void nec_ground::setup_cliff(nec_float in_eprs2,
      nec_float in_sig2,
      nec_float clt, nec_float cht)
{
  cliff_edge_distance = clt;
  cliff_height = cht;
  epsr2 = in_eprs2;
  sig2 = in_sig2;
}


#include <ostream>
using namespace std;

#include "electromag.h"
#include "nec_results.h"

nec_complex nec_ground::get_zrati2(nec_float _wavelength)  {
  return sqrt(1.0 / nec_complex(epsr2,-sig2 * _wavelength * em::impedance_over_2pi()));
}


void nec_ground::write_to_file_aux(std::ostream& os, output_helper& oh, int _ifar)
{
  if ( _ifar >= 2)  {
    oh.section_start("FAR FIELD GROUND PARAMETERS");
    
    if ( _ifar > 3) {
      os << endl;
      os << "                                        RADIAL WIRE GROUND SCREEN" << endl;
      os << "                                        "; oh.int_out(5, radial_wire_count); os << " WIRES" << endl;
      os << "                                        WIRE LENGTH= "; oh.real_out(8,2, radial_wire_length,false); os << " METERS" << endl;
      os << "                                        WIRE RADIUS= "; oh.real_out(10,3, radial_wire_radius); os << " METERS" << endl;
    } /* if ( _ifar > 3) */
    
    if ( _ifar != 4 )  {
      std::string hclif;
      if ( (_ifar == 2) || (_ifar == 5) )
        hclif = "LINEAR";
      if ( (_ifar == 3) || (_ifar == 6) )
        hclif= "CIRCLE";

      os << endl;
      os << "                                        " << hclif << " CLIFF" << endl;
      os << "                                        EDGE DISTANCE= "; oh.real_out(9,2,cliff_edge_distance,false); os << " METERS" << endl;
      os << "                                        HEIGHT= "; oh.real_out(8,2,cliff_height,false); os << " METERS" << endl;
      os << "                                        SECOND MEDIUM -" << endl;
      os << "                                        RELATIVE DIELECTRIC CONST.= "; oh.real_out(7,3,epsr2, false); os << endl;
      os << "                                        CONDUCTIVITY= "; oh.real_out(10,3,sig2,false); os << " MHOS" << endl;
    } /* if ( _ifar != 4 ) */
  } /* if ( _ifar >= 2) */
}

void nec_ground::calculate_antenna_environment(c_ground_wave& ground_wave, nec_float freq_mhz)
{
  if ( false == present()) {
      return;
  }
  
  frati=cplx_10();

  if (type_perfect()) {
    return;
  }
  
  nec_float wavelength = em::get_wavelength(freq_mhz * 1.0e6);
  
  if ( sig < 0.)
    sig = - sig/(em::impedance_over_2pi()*wavelength);

  _epsc = nec_complex( epsr, -sig*wavelength*em::impedance_over_2pi());
  zrati = 1.0/ sqrt( _epsc);
  
  ground_wave.set_u(zrati);

  if (  radial_wire_count != 0) {
    scrwl=  radial_wire_length/ wavelength;
    scrwr=  radial_wire_radius/ wavelength;
    m_t1 = cplx_01()*2367.067/ (nec_float) radial_wire_count;
    t2 = scrwr * (nec_float) radial_wire_count;
  }

  if (type_sommerfeld_norton()) {
    // calculate the Sommerfeld Norton ground stuff.
    _ggrid.sommerfeld( epsr, sig, wavelength );
    
    frati = (_epsc-1.0)/(_epsc+1.0);
    if ( abs(( _ggrid.m_epscf- _epsc)/ _epsc) >= 1.0e-3 ) {
      nec_stop("ERROR IN GROUND PARAMETERS -"
        "\n COMPLEX DIELECTRIC CONSTANT FROM FILE IS: %12.5E%+12.5Ej"
        "\n REQUESTED: %12.5E%+12.5Ej",
        real(_ggrid.m_epscf), imag(_ggrid.m_epscf), 
        real(_epsc), imag(_epsc) );
    }
  }
}

void nec_ground::output_antenna_environment(nec_output_file& _output)
{
  _output.end_section();
  _output.line("                            -------- ANTENNA ENVIRONMENT --------" );
  
  if ( false == present()) {
    _output.line("                            FREE SPACE" );
    return;
  }
  
  if (type_perfect()) {
    _output.line("                            PERFECT GROUND" );
    return;
  }
  
  if (  radial_wire_count != 0) {
    _output.line(
        "                            RADIAL WIRE GROUND SCREEN");
    _output.nec_printf(
        "                            %d WIRES\n"
        "                            WIRE LENGTH: %8.2f METERS\n"
        "                            WIRE RADIUS: %10.3E METERS",
        radial_wire_count,  radial_wire_length,  radial_wire_radius );

    _output.endl();
    _output.line("                            MEDIUM UNDER SCREEN -" );
  }

  if (false == type_sommerfeld_norton()) {
    _output.line("                            FINITE GROUND - REFLECTION COEFFICIENT APPROXIMATION" );
  } else {
    _output.line("                            FINITE GROUND - SOMMERFELD SOLUTION" );
  }

  _output.endl();
  _output.nec_printf(
      "                            "
      "RELATIVE DIELECTRIC CONST: %.3f\n"
      "                            "
      "CONDUCTIVITY: %10.3E MHOS/METER\n"
      "                            "
      "COMPLEX DIELECTRIC CONSTANT: %11.4E%+11.4Ej",
      epsr, sig, real(_epsc), imag(_epsc) );
}
