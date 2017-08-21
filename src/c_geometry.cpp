/*
  Copyright (C) 2004-2011  Timothy C.A. Molteno
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
#include "c_geometry.h"

#include "nec_context.h"
#include "nec_exception.h"

#include <cstring>
#include <stdint.h>

c_geometry::c_geometry() 
    : patch_x1(0,0,0), patch_x2(0,0,0), patch_x3(0,0,0), patch_x4(0,0,0)  {
  n_segments = 0;
  np = 0;   // n_segments is the number of segments
  
  m = 0;
  mp = 0;    // m is the number of patches
  
  m_ipsym = 0;
  
  n_plus_2m = 0;
  n_plus_3m = 0;
  
  jsno = 0;
  nscon = 0;
  maxcon = 0;
  
  m_context = NULL;
  m_output = NULL;
}

void c_geometry::set_context(nec_context* in_context)  {
  m_context = in_context;
  m_output = &m_context->m_output;
}

/*! \brief Get a segment number for a specified tag.
  \param in_tag The tag
  \param in_m The mth segment with the specified tag will be returned.
  \return The segment number of the mth segment having the
  tag number in_tag.  if in_tag=0 segment number m is returned.
*/
int c_geometry::get_segment_number( int in_tag, int in_m)
{
  ASSERT(in_tag >= 0);
  ASSERT(in_m >= 0);
  
  
  if (in_m <= 0)
  {
    throw new nec_exception("CHECK DATA, PARAMETER SPECIFYING SEGMENT POSITION IN A GROUP OF EQUAL TAGS MUST NOT BE ZERO" );
  }
  
  if ( 0 == in_tag)
  {
    return( in_m );
  }
  
  int tag_seg_count=0;
  for (int i = 0; i < n_segments; i++ )
  {
    if ( segment_tags[i] == in_tag )
    {
      tag_seg_count++;
      if ( tag_seg_count == in_m)
      {
        return( i+1 );
      }
    }
  }
  
  throw new nec_exception("NO SEGMENT HAS AN ITAG OF ", in_tag);
  
  return 0;
}

#include "c_plot_card.h"

#include <algorithm>

void str_toupper(std::string &str);
void str_toupper(std::string &str)
{
  std::transform(str.begin(), 
    str.end(), 
    str.begin(),
    ::toupper);
}
void c_geometry::parse_geometry(nec_context* in_context, FILE* input_fp )
{
  char gm[3];
  const char ipt[4] = { 'P', 'R', 'T', 'Q' };
  
  /* input card mnemonic list */
  /* "XT" stands for "exit", added for testing */
/*  #define GM_NUM  12
  char *atst[GM_NUM] =
  {
    "GW", "GX", "GR", "GS", "GE", "GM", \
    "SP", "SM", "GA", "SC", "GH", "GF"
  };*/
  
  bool print_structure_spec = true;
  
  int nwire, isct, i1, i2;
  int card_int_1, card_int_2; /* The two integer parameters from the geometry card */
  nec_float rad, xs1, xs2, ys1, ys2, zs1, zs2, x4=0, y4=0, z4=0;
  nec_float x3=0, y3=0, z3=0, xw1, xw2, yw1, yw2, zw1, zw2;
  nec_float dummy;
  
  m_ipsym=0;
  nwire=0;
  n_segments=0;
  np=0;
  m=0;
  mp=0;
  isct=0;
  
  /* read geometry data card and branch to */
  /* section for operation requested */
  do
  {
    read_geometry_card(input_fp, gm, &card_int_1, &card_int_2, &xw1, &yw1, &zw1, &xw2, &yw2, &zw2, &rad);
  
    /* identify card id mnemonic */
    std::string card_id(gm);
    str_toupper(card_id);
//    int gm_num;
//    for( gm_num = 0; gm_num < GM_NUM; gm_num++ )
//      if ( strncmp( gm, atst[gm_num], 2) == 0 )
//        break;

    if ( print_structure_spec )
    {
      m_output->end_section();
      m_output->set_indent(32);
      m_output->line("-------- STRUCTURE SPECIFICATION --------");
      m_output->line("COORDINATES MUST BE INPUT IN" );
      m_output->line("METERS OR BE SCALED TO METERS" );
      m_output->line("BEFORE STRUCTURE INPUT IS ENDED" );
      m_output->set_indent(0);
    
      m_output->line("  WIRE                                                                                 SEG FIRST  LAST  TAG");
      m_output->line("   No:        X1         Y1         Z1         X2         Y2         Z2       RADIUS   No:   SEG   SEG  No:" );
    
      print_structure_spec = false;
    }

    if ( card_id != "SC") // gm_num != 10 )
      isct=0;

    /* "gw" card, generate segment data for straight wire.
      GW    STRAIGHT WIRE, ENDS 1,2
        card_int_1- TAG NO.
        card_int_2- NO. SEGMENTS
        xw1- X1
        F2- Y1
        F3- Z1
        F4- X2
        F5- Y2
        F6- Z2
        F7- WIRE RAD., 0=USE GC FOR TAPERED WIRE
    */
    if (card_id == "GW")
    {
      int wire_segment_count = card_int_2;
      int wire_tag = card_int_1;
      
      nwire++;
    
      // output some wire diagnostics.
      m_output->nec_printf( "\n"
        " %5d  %10.4f %10.4f %10.4f %10.4f"
        " %10.4f %10.4f %10.4f %5d %5d %5d %4d",
        nwire, xw1, yw1, zw1, xw2, yw2, zw2, rad, wire_segment_count, n_segments+1, n_segments + wire_segment_count, wire_tag );
    
      if ( rad != 0.0)  // rad == 0 implies a tapered wire
      {
        xs1 = 1.0;
        ys1 = 1.0;
      }
      else
      {
        int ix,iy;
        read_geometry_card(input_fp, gm, &ix, &iy, &xs1, &ys1, &zs1,
          &dummy, &dummy, &dummy, &dummy);
      
        if ( strcmp(gm, "GC" ) != 0 )
        {
          throw new nec_exception("GEOMETRY DATA CARD ERROR" );
        }
      
        m_output->nec_printf(
          "\n  ABOVE WIRE IS TAPERED.  SEGMENT LENGTH RATIO: %9.5f\n"
          "                                 "
          "RADIUS FROM: %9.5f TO: %9.5f", xs1, ys1, zs1 );
      
        if ( (ys1 == 0) || (zs1 == 0) )
        {
          throw new nec_exception("GEOMETRY DATA CARD ERROR" );
        }
      
        rad= ys1;
        ys1= pow( (zs1/ys1), (1./(wire_segment_count-1.)) );
      }
    
      wire(wire_tag, wire_segment_count, xw1, yw1, zw1, xw2, yw2, zw2, rad, xs1, ys1);
    }

    /* reflect structure along x,y, or z */
    /* axes or rotate to form cylinder.  */
    /* "gx" card */
    else if (card_id == "GX")
    {  gx_card(card_int_1, card_int_2);
    }
  
    /* "gr" card */
    else if (card_id == "GR")
    {
    
      m_output->nec_printf(
        "\n  STRUCTURE ROTATED ABOUT Z-AXIS %d TIMES"
        " - LABELS INCREMENTED BY %d\n", card_int_2, card_int_1 );
    
      int ix = -1;
      int iy = 0;
      int iz = 0;
      reflect( ix, iy, iz, card_int_1, card_int_2);
    }
  
    /* "gs" card, scale structure dimensions by factor xw1. */
    else if (card_id == "GS")
    {      
      m_output->nec_printf(
        "\n     STRUCTURE SCALED BY FACTOR: %10.5f", xw1 );
      
      scale(xw1);  
    }

    /* "ge" card, terminate structure geometry input. */
    else if (card_id == "GE")
    {
      geometry_complete(in_context, card_int_1);
      return;
    }

    /* "gm" card, move structure or reproduce */
    /* original structure in new positions.   */
    else if (card_id == "GM")
    {
      m_output->nec_printf(
        "\n     THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS:\n"
        "   %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
        card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, rad );
    
      xw1= degrees_to_rad(xw1);
      yw1= degrees_to_rad(yw1);
      zw1= degrees_to_rad(zw1);
    
      move( xw1, yw1, zw1, xw2, yw2, zw2, (int)( rad+.5), card_int_2, card_int_1);
    }
    
    /* "sp" card, generate single new patch */
    else if (card_id == "SP")
    {
      i1= m+1;
      card_int_2++;
    
      if ( card_int_1 != 0)
      {
        throw new nec_exception("PATCH DATA ERROR" );
      }
    
      m_output->nec_printf( "\n"
        " %5d%c %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
        i1, ipt[card_int_2-1], xw1, yw1, zw1, xw2, yw2, zw2 );
    
      if ( (card_int_2 == 2) || (card_int_2 == 4) )
        isct=1;
    
      if ( card_int_2 > 1)
      {
        // read another geometry card for the rest of the patch data. 
        int ix,iy;
        read_geometry_card(input_fp, gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);
      
        if ( (card_int_2 == 2) || (card_int_1 > 0) )
        {
          x4= xw1+ x3- xw2;
          y4= yw1+ y3- yw2;
          z4= zw1+ z3- zw2;
        }
      
        m_output->nec_printf( "\n"
          "      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
          x3, y3, z3, x4, y4, z4 );
      
        if ( strcmp(gm, "SC") != 0 )
        {
          throw new nec_exception("PATCH DATA ERROR" );
        }
      } /* if ( card_int_2 > 1) */
      else
      {
        xw2= degrees_to_rad(xw2);
        yw2= degrees_to_rad(yw2);
      }
    
      patch( card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);
    }
  
    /* "sm" card, generate multiple-patch surface */
    else if (card_id == "SM")
    {
      i1= m+1;
      m_output->nec_printf( "\n"
        " %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f"
        "     SURFACE - %d BY %d PATCHES",
        i1, ipt[1], xw1, yw1, zw1, xw2, yw2, zw2, card_int_1, card_int_2 );
    
      if ( (card_int_1 < 1) || (card_int_2 < 1) )
      {
        throw new nec_exception("PATCH DATA ERROR" );
      }
    
      int ix,iy;
      read_geometry_card(input_fp, gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);
    
      if ( (card_int_2 == 2) || (card_int_1 > 0) )
      {
        x4= xw1+ x3- xw2;
        y4= yw1+ y3- yw2;
        z4= zw1+ z3- zw2;
      }
    
      m_output->nec_printf( "\n"
        "      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
        x3, y3, z3, x4, y4, z4 );
    
      if ( strcmp(gm, "SC" ) != 0 )
      {
        throw new nec_exception("PATCH DATA ERROR" );
      }
    
      patch( card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);
    }
    
    /* "ga" card, generate segment data for wire arc */
    else if (card_id == "GA")
    {
      nwire++;
      i1= n_segments+1;
      i2= n_segments+ card_int_2;
    
      m_output->nec_printf( "\n"
        " %5d  ARC RADIUS: %9.5f  FROM: %8.3f TO: %8.3f DEGREES"
        "       %11.5f %5d %5d %5d %4d",
        nwire, xw1, yw1, zw1, xw2, card_int_2, i1, i2, card_int_1 );
    
      arc( card_int_1, card_int_2, xw1, yw1, zw1, xw2);
    }
  
    /* "sc" card */
    else if (card_id == "SC")
    {
      if ( isct == 0)
      {
        throw new nec_exception("PATCH DATA ERROR" );
      }
    
      i1= m+1;
      card_int_2++;
    
      if ( (card_int_1 != 0) || ((card_int_2 != 2) && (card_int_2 != 4)) )
      {
        throw new nec_exception("PATCH DATA ERROR" );
      }
    
      xs1= x4;
      ys1= y4;
      zs1= z4;
      xs2= x3;
      ys2= y3;
      zs2= z3;
      x3= xw1;
      y3= yw1;
      z3= zw1;
    
      if ( card_int_2 == 4)
      {
        x4= xw2;
        y4= yw2;
        z4= zw2;
      }
    
      xw1= xs1;
      yw1= ys1;
      zw1= zs1;
      xw2= xs2;
      yw2= ys2;
      zw2= zs2;
    
      if ( card_int_2 != 4)
      {
        x4= xw1+ x3- xw2;
        y4= yw1+ y3- yw2;
        z4= zw1+ z3- zw2;
      }
    
      m_output->nec_printf( "\n"
        " %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
        i1, ipt[card_int_2-1], xw1, yw1, zw1, xw2, yw2, zw2 );
    
      m_output->nec_printf( "\n"
        "      %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f",
        x3, y3, z3, x4, y4, z4 );
    
      patch( card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);
    }
  
    /* "gh" card, generate helix */
    else if (card_id == "GH")
    {
      nwire++;
      i1= n_segments+1;
      i2= n_segments+ card_int_2;
    
      m_output->nec_printf( "\n"
        " %5d HELIX STRUCTURE - SPACING OF TURNS: %8.3f AXIAL"
        " LENGTH: %8.3f  %8.3f %5d %5d %5d %4d\n      "
        " RADIUS X1:%8.3f Y1:%8.3f X2:%8.3f Y2:%8.3f ",
        nwire, xw1, yw1, rad, card_int_2, i1, i2, card_int_1, zw1, xw2, yw2, zw2 );
      int tag_id(card_int_1);
      int segment_count(card_int_2);
      nec_float s(xw1);
      nec_float hl(yw1);
      nec_float a1(zw1);
      nec_float b1(xw2);
      nec_float a2(yw2);
      nec_float b2(zw2);
      
      helix(tag_id, segment_count,
            s, hl, a1, b1, a2, b2, rad);
    }

    /* "gf" card, not supported */
    else if (card_id == "GF")
      throw new nec_exception("NGF solution option not supported");
  
    /* error message */
    else
    {
      m_output->nec_printf( "\n  GEOMETRY DATA CARD ERROR" );
      m_output->nec_printf( "\n"
        " %2s %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
        gm, card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, rad );
    
      throw new nec_exception("GEOMETRY DATA CARD ERROR");
    }
 
  } /* do */
  while( true );
}

#include "nec_wire.h"
/**
  We have finished with the geometry description, now connect 
  things up.
*/
void c_geometry::geometry_complete(nec_context* in_context, int gpflag)
{
  if (0 == np + mp)
    throw new nec_exception("Geometry has no wires or patches.");
    
  /* Check to see whether any wires intersect with one another */
  for (uint32_t i=0; i<m_wires.size(); i++)
  {
    nec_wire a = m_wires[i];
    for (uint32_t j=0; j<m_wires.size(); j++)
    {
      if (i > j)
      {
        nec_wire b = m_wires[j];
        vector<nec_wire> wires = a.intersect(b);
        if (wires.size() > 2)
        {
          nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR -- WIRE #");
          nex->append(j+1);
          nex->append(" (TAG ID #"); nex->append(b.tag_id());
          nex->append(") INTERSECTS WIRE #");
          nex->append(i+1);
          nex->append(" (TAG ID #"); nex->append(a.tag_id()); nex->append(")");
          throw nex;
        }
      }
    }
  }

// now proceed and complete the geometry setup...
  // Check here that patches form a closed surfaceAntennaInput
  /*
    From Jerry Burke (the original author of NEC2):
    Patches are modeled in NEC-2 
    and NEC-4 with the magnetic field integral equation, and it is valid 
    only for closed perfectly conducting surfaces. �You start with a 
    single patch and connect a wire to its center so that it splits into 
    four patches. �You would need to use the SM command to make a surface 
    with more patches, then use other SM or MV commands to form the other 
    five faces of a closed box. �People often do 
    use the NEC patch model for non-closed surfaces, and it is completely 
    invalid. �It would be very useful to have a check for non-closed 
    patch surfaces in NEC-2 or in GUIs written for it , but determining 
    whether a bunch of patches forms a closed surface is not easy. �If 
    you want to model a monopole on a surface in NEC-2, you should model 
    the surface as a wire grid.
  */

   connect_segments( gpflag);

  if ( n_segments != 0)
  {
    /* Allocate wire buffers */
    segment_length.resize(n_segments);
    sab.resize(n_segments);
    cab.resize(n_segments);
    salp.resize(n_segments);
  
    m_output->nec_printf( "\n\n\n"
      "                              "
      " ---------- SEGMENTATION DATA ----------\n"
      "                                       "
      " COORDINATES IN METERS\n"
      "                           "
      " I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I\n" );
  
    m_output->nec_printf( "\n"
      "   SEG    COORDINATES OF SEGM CENTER     SEGM    ORIENTATION"
      " ANGLES    WIRE    CONNECTION DATA   TAG\n"
      "   No:       X         Y         Z      LENGTH     ALPHA     "
      " BETA    RADIUS    I-     I    I+   NO:" );
  
    for(int i = 0; i < n_segments; i++ )
    {
      nec_float xw1= x2[i]- x[i];
      nec_float yw1= y2[i]- y[i];
      nec_float zw1= z2[i]- z[i];
      x[i]=( x[i]+ x2[i])*.5;
      y[i]=( y[i]+ y2[i])*.5;
      z[i]=( z[i]+ z2[i])*.5;
      
      nec_float xw2= xw1* xw1+ yw1* yw1+ zw1* zw1;
      nec_float yw2= sqrt( xw2);
      yw2=( xw2/ yw2+ yw2)*.5;
      segment_length[i]= yw2;
      cab[i]= xw1/ yw2;
      sab[i]= yw1/ yw2;
      xw2= zw1/ yw2;
  
      if ( xw2 > 1.)
        xw2=1.;
      if ( xw2 < -1.)
        xw2 = -1.;
  
      salp[i]= xw2;
      xw2= rad_to_degrees(asin( xw2));
      yw2= rad_to_degrees(atan2( yw1, xw1));
  
      m_output->nec_printf( "\n"
        " %5d %9.4f %9.4f %9.4f %9.4f"
        " %9.4f %9.4f %9.4f %5d %5d %5d %5d",
        i+1, x[i], y[i], z[i], segment_length[i], xw2, yw2,
        segment_radius[i], icon1[i], i+1, icon2[i], segment_tags[i] );
  
      in_context->plot_card.plot_segments(i,x,y,z,segment_length,xw2,yw2,segment_radius,icon1,icon2);
  
      if ( (segment_length[i] <= 1.e-20) || (segment_radius[i] <= 0.) )
      {
        throw new nec_exception("SEGMENT DATA ERROR" );
      }
    } /* for( i = 0; i < n_segments; i++ ) */
  } /* if ( n_segments != 0) */

  if ( m != 0)
  {
    m_output->nec_printf( "\n\n\n"
      "                                   "
      " --------- SURFACE PATCH DATA ---------\n"
      "                                            "
      " COORDINATES IN METERS\n\n"
      " PATCH      COORD. OF PATCH CENTER           UNIT NORMAL VECTOR      "
      " PATCH           COMPONENTS OF UNIT TANGENT VECTORS\n"
      "  NO:       X          Y          Z          X        Y        Z      "
      " AREA         X1       Y1       Z1        X2       Y2      Z2" );
  
    for(int i = 0; i < m; i++ )
    {
      nec_float xw1=( t1y[i]* t2z[i]- t1z[i]* t2y[i])* psalp[i];
      nec_float yw1=( t1z[i]* t2x[i]- t1x[i]* t2z[i])* psalp[i];
      nec_float zw1=( t1x[i]* t2y[i]- t1y[i]* t2x[i])* psalp[i];
  
      m_output->nec_printf( "\n"
        " %4d %10.5f %10.5f %10.5f  %8.4f %8.4f %8.4f"
        " %10.5f  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f",
        i+1, px[i], py[i], pz[i], xw1, yw1, zw1, pbi[i],
        t1x[i], t1y[i], t1z[i], t2x[i], t2y[i], t2z[i] );
    } /* for( i = 0; i < m; i++ ) */
  } /* if ( m == 0) */

  n_plus_2m = n_segments+2*m;
  n_plus_3m = n_segments+3*m;

  x_unscaled.resize(n_segments);
  y_unscaled.resize(n_segments);
  z_unscaled.resize(n_segments);
  si_unscaled.resize(n_segments);
  bi_unscaled.resize(n_segments);
  
  px_unscaled.resize(m);
  py_unscaled.resize(m);
  pz_unscaled.resize(m);
  pbi_unscaled.resize(m);
  
  // Fill the unscaled segments...
  for (int i = 0; i < n_segments; i++ )
  {
    x_unscaled[i]= x[i];
    y_unscaled[i]= y[i];
    z_unscaled[i]= z[i];
    si_unscaled[i]= segment_length[i];
    bi_unscaled[i]= segment_radius[i];
  }
  // Fill the unscaled patches...
  for (int i = 0; i < m; i++ )
  {
    px_unscaled[i]= px[i];
    py_unscaled[i]= py[i];
    pz_unscaled[i]= pz[i];
    pbi_unscaled[i]= pbi[i];
  }
}


/*! \brief Generates segment geometry for a straingt wire
  \param tag_id
  \param segment_count Number of Elements (should be around 12-20 per wavelength)
  \param rad Wire radius of first segment (in Meters)
  \param rdel Ratio of the length of a segment to the length of the previous segment.  (Set to 1.0 if segments have uniform length)
  \param rrad The ratio of the radii of adjacent segments (Set to 1.0 if not tapered)
*/
void c_geometry::wire( int tag_id, int segment_count, nec_float xw1, nec_float yw1, nec_float zw1,
    nec_float xw2, nec_float yw2, nec_float zw2, nec_float rad,
    nec_float rdel, nec_float rrad )
{
  nec_float delz, rd, fns, radz;
  
  int istart = n_segments;
  n_segments += segment_count;
  np= n_segments;
  mp= m;
  m_ipsym=0;
  
  if ( segment_count < 1)
    return;
  
  /* Reallocate tags buffer */
  segment_tags.resize(n_segments + m);
  
  /* Reallocate wire buffers */
  x.resize(n_segments);
  y.resize(n_segments);
  z.resize(n_segments);
  x2.resize(n_segments);
  y2.resize(n_segments);
  z2.resize(n_segments);
  segment_radius.resize(n_segments);
  
  nec_3vector dx(xw2- xw1, yw2- yw1, zw2- zw1);
  
  if ( fabs(rdel-1.0) >= 1.0e-6) // Use a tapered wire
  {
    delz= dx.norm();
    dx /= delz;
    delz= delz*(1.- rdel)/(1.- pow(rdel, segment_count) );
    rd= rdel;
  }
  else
  {
    fns= segment_count;
    dx /= fns;
    delz=1.0;
    rd=1.0;
  }
  
  /*
  There is no restriction on the angle between two wires, but accuracy will be lost if the center of a segment falls within the
  volume of the wire the segment connects to. The risk of this reduces as the angle between wires approaches 180 degrees.
  
  Wires which intersect away from their ends are not connected, but errors will occur if one wire occupies the space of another one. For accuracy, separate wire centers by several radii of the largest wire. 
  */
  // check that none of the existing wires intersect with the midpoint
  // of the first and last segment
  nec_3vector wire_start(xw1,yw1,zw1);
  nec_3vector wire_end(xw2,yw2,zw2);
  
  nec_3vector seg_midpoint(wire_start + (dx/2)*delz);
  nec_3vector end_seg_midpoint = wire_end - (dx*delz / 2);
  /* Check to see whether any wires intersect with the segment_midpoint */
  for (uint32_t i=0; i<m_wires.size(); i++)
  {
    nec_wire a = m_wires[i];
    if (a.intersect(seg_midpoint))
    {
      nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR -- FIRST SEGMENT MIDPOINT");
      nex->append(" OF WIRE #");
      nex->append(m_wires.size()+1);
      nex->append(" (TAG ID #"); nex->append(tag_id);
      nex->append(") INTERSECTS WIRE #");
      nex->append(i+1);
      nex->append(" (TAG ID #"); nex->append(a.tag_id()); nex->append(")");
      
      throw nex;
    }
    if (a.intersect(end_seg_midpoint))
    {
      nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR -- LAST SEGMENT MIDPOINT");
      nex->append(" OF WIRE #");
      nex->append(m_wires.size()+1);
      nex->append(" (TAG ID #"); nex->append(tag_id);
      nex->append(") INTERSECTS WIRE #");
      nex->append(i+1);
      nex->append(" (TAG ID #"); nex->append(a.tag_id()); nex->append(")");
      throw nex;
    }
  }
  

  radz= rad;
  nec_3vector xs1(xw1,yw1,zw1);
  nec_3vector x_end(xw2,yw2,zw2);

  m_wires.push_back(nec_wire(xs1, x_end, rad, tag_id));

  for (int i = istart; i < n_segments; i++ )
  {
    segment_tags[i]= tag_id;
    nec_3vector xs2(xs1 + dx*delz);
    x[i]= xs1.x();
    y[i]= xs1.y();
    z[i]= xs1.z();
    x2[i]= xs2.x();
    y2[i]= xs2.y();
    z2[i]= xs2.z();
    ASSERT(0.0 != radz);
    segment_radius[i]= radz;
    delz= delz* rd;
    radz= radz* rrad;
    xs1 = xs2;
  }


  x2[n_segments-1]= xw2;
  y2[n_segments-1]= yw2;
  z2[n_segments-1]= zw2;

}

/*-----------------------------------------------------------------------*/

/* subroutine helix generates segment geometry
  data for a helix of segment_count segments 
 
           S   (F1) - Spacing between turns.
           HL  (F2) - Total length of the helix.
           A1  (F3) - Radius in x at z = 0.
           B1  (F4) - Radius in y at z = 0.
           A2  (F5) - Radius in x at z = HL.
           B2  (F6) - Radius in y at z = HL.
           RAD (F7) - Radius of wire.
*/
void c_geometry::helix(int tag_id, int segment_count, nec_float s, nec_float hl, nec_float a1, nec_float b1,
    nec_float a2, nec_float b2, nec_float rad)
{
  int ist;
  nec_float zinc, sangle, hdia, turn, pitch, hmaj, hmin;
  
  ist= n_segments;
  n_segments += segment_count;
  np= n_segments;
  mp= m;
  m_ipsym=0;
  
  if ( segment_count < 1)
    return;
  
  zinc= fabs( hl/ segment_count);
  
  segment_tags.resize(n_segments+m); /*????*/
  
  /* Reallocate wire buffers */
  x.resize(n_segments);
  y.resize(n_segments);
  z.resize(n_segments);
  x2.resize(n_segments);
  y2.resize(n_segments);
  z2.resize(n_segments);
  segment_radius.resize(n_segments);
  
  z[ist]=0.;
  for(int i = ist; i < n_segments; i++ ) {
    segment_radius[i]= rad;
    segment_tags[i]= tag_id;
  
    if ( i != ist )
    z[i]= z[i-1]+ zinc;
  
    z2[i]= z[i]+ zinc;
  
    if ( a2 == a1) {
      if ( b1 == 0.)
        b1= a1;
    
      x[i]= a1* cos(2.* pi()* z[i]/ s);
      y[i]= b1* sin(2.* pi()* z[i]/ s);
      x2[i]= a1* cos(2.* pi()* z2[i]/ s);
      y2[i]= b1* sin(2.* pi()* z2[i]/ s);
    } else {
      if ( b2 == 0.)
        b2= a2;
    
      x[i]=( a1+( a2- a1)* z[i]/ fabs( hl))* cos(2.* pi()* z[i]/ s);
      y[i]=( b1+( b2- b1)* z[i]/ fabs( hl))* sin(2.* pi()* z[i]/ s);
      x2[i]=( a1+( a2- a1)* z2[i]/ fabs( hl))* cos(2.* pi()* z2[i]/ s);
      y2[i]=( b1+( b2- b1)* z2[i]/ fabs( hl))* sin(2.* pi()* z2[i]/ s);  
    } /* if ( a2 == a1) */
  
    if ( hl <= 0.) {  
      nec_float copy= x[i];
      x[i]= y[i];
      y[i]= copy;
      copy= x2[i];
      x2[i]= y2[i];
      y2[i]= copy;
    }
  } /* for( i = ist; i < n_segments; i++ ) */
  
  if ( a2 != a1) {
    sangle= atan( a2/( fabs( hl)+( fabs( hl)* a1)/( a2- a1)));
    m_output->nec_printf(
      "\n       THE CONE ANGLE OF THE SPIRAL IS %10.4f", sangle );
    return;
  }
  
  if ( a1 == b1) {
    hdia=2.* a1;
    turn= hdia* pi();
    pitch= atan( s/( pi()* hdia));
    turn= turn/ cos( pitch);
    pitch=180.* pitch/ pi();
  } else {
    if ( a1 >= b1) {
      hmaj=2.* a1;
      hmin=2.* b1;
    } else {
      hmaj=2.* b1;
      hmin=2.* a1;
    }
  
    hdia= sqrt(( hmaj*hmaj+ hmin*hmin)/2* hmaj);
    turn=2.* pi()* hdia;
    pitch=(180./ pi())* atan( s/( pi()* hdia));
  } /* if ( a1 == b1) */
  
  m_output->nec_printf( "\n"
    "       THE PITCH ANGLE IS: %.4f    THE LENGTH OF WIRE/TURN IS: %.4f",
    pitch, turn );
}

/*-----------------------------------------------------------------------*/


/* subroutine move moves the structure with respect to its */
/* coordinate system or reproduces structure in new positions. */
/* structure is rotated about x,y,z axes by rox,roy,roz */
/* respectively, then shifted by xs,ys,zs */
void c_geometry::move( nec_float rox, nec_float roy, nec_float roz, nec_float xs,
    nec_float ys, nec_float zs, int its, int nrpt, int itgi )
{
  DEBUG_TRACE("move " << nrpt << " Copies");
  int nrp, ix, i1, k;
  nec_float sps, cps, sth, cth, sph, cph, xx, xy;
  nec_float xz, yx, yy, yz, zx, zy, zz, xi, yi, zi;

  if ( fabs( rox)+ fabs( roy) > 1.0e-10)
    m_ipsym= m_ipsym*3;

  sps= sin( rox);
  cps= cos( rox);
  sth= sin( roy);
  cth= cos( roy);
  sph= sin( roz);
  cph= cos( roz);
  xx= cph* cth;
  xy= cph* sth* sps- sph* cps;
  xz= cph* sth* cps+ sph* sps;
  yx= sph* cth;
  yy= sph* sth* sps+ cph* cps;
  yz= sph* sth* cps- cph* sps;
  zx = - sth;
  zy= cth* sps;
  zz= cth* cps;

  if ( nrpt == 0)
    nrp=1;
  else
    nrp= nrpt;

  ix=1;
  if ( n_segments > 0) {
    i1= get_segment_number( its, 1);
    if ( i1 < 1)
      i1= 1;

    ix= i1;
    if ( nrpt == 0) {
      k= i1-1;
    } else {
      k= n_segments;
      /* Reallocate tags buffer */
      segment_tags.resize(n_segments+m + (n_segments+1-i1)*nrpt);
      // mreq = n_segments+m + (n_segments+1-i1)*nrpt;
      // segment_tags.resize(mreq);

      /* Reallocate wire buffers */
      int new_size = (n_segments+(n_segments+1-i1)*nrpt);
      x.resize(new_size);
      y.resize(new_size);
      z.resize(new_size);
      x2.resize(new_size);
      y2.resize(new_size);
      z2.resize(new_size);
      segment_radius.resize(new_size);
    }

    for (int ir = 0; ir < nrp; ir++ ) {
      DEBUG_TRACE("GM: Segment Copy #" << ir);
      for (int i = i1-1; i < n_segments; i++ ) {
        xi= x[i];
        yi= y[i];
        zi= z[i];
        x[k]= xi* xx+ yi* xy+ zi* xz+ xs;
        y[k]= xi* yx+ yi* yy+ zi* yz+ ys;
        z[k]= xi* zx+ yi* zy+ zi* zz+ zs;
        xi= x2[i];
        yi= y2[i];
        zi= z2[i];
        x2[k]= xi* xx+ yi* xy+ zi* xz+ xs;
        y2[k]= xi* yx+ yi* yy+ zi* yz+ ys;
        z2[k]= xi* zx+ yi* zy+ zi* zz+ zs;
        segment_radius[k]= segment_radius[i];
        segment_tags[k]= segment_tags[i];
        if ( segment_tags[i] != 0)
          segment_tags[k]= segment_tags[i]+ itgi;

        k++;
      } /* for( i = i1; i < n_segments; i++ ) */

      i1= n_segments+1;
      n_segments= k;
    } /* for( ir = 0; ir < nrp; ir++ ) */
  } /* if ( n_segments >= n2) */

  if ( m > 0) {
    i1 = 0;
    if ( nrpt == 0)
      k= 0;
    else
      k = m;

    /* Reallocate patch buffers */
    int new_size = m * (1+nrpt);
    px.resize(new_size);
    py.resize(new_size);
    pz.resize(new_size);
    t1x.resize(new_size);
    t1y.resize(new_size);
    t1z.resize(new_size);
    t2x.resize(new_size);
    t2y.resize(new_size);
    t2z.resize(new_size);
    pbi.resize(new_size);
    psalp.resize(new_size);

    for (int ii = 0; ii < nrp; ii++ ) {
      DEBUG_TRACE("GM: Patch Copy #" << ii);
      for(int i = i1; i < m; i++ ) {
        xi= px[i];
        yi= py[i];
        zi= pz[i];
        px[k]= xi* xx+ yi* xy+ zi* xz+ xs;
        py[k]= xi* yx+ yi* yy+ zi* yz+ ys;
        pz[k]= xi* zx+ yi* zy+ zi* zz+ zs;
        xi= t1x[i];
        yi= t1y[i];
        zi= t1z[i];
        t1x[k]= xi* xx+ yi* xy+ zi* xz;
        t1y[k]= xi* yx+ yi* yy+ zi* yz;
        t1z[k]= xi* zx+ yi* zy+ zi* zz;
        xi= t2x[i];
        yi= t2y[i];
        zi= t2z[i];
        t2x[k]= xi* xx+ yi* xy+ zi* xz;
        t2y[k]= xi* yx+ yi* yy+ zi* yz;
        t2z[k]= xi* zx+ yi* zy+ zi* zz;
        psalp[k]= psalp[i];
        pbi[k]= pbi[i];
        k++;
      } /* for( i = i1; i < m; i++ ) */
      i1= m;
      m = k;
    } /* for( ii = 0; ii < nrp; ii++ ) */
  } /* if ( m >= m2) */

  if ( (nrpt == 0) && (ix == 1) )
    return;

  np= n_segments;
  mp= m;
  m_ipsym=0;
}

/*-----------------------------------------------------------------------*/

/* reflects partial structure along x,y, or z axes or rotates */
/* structure to complete a symmetric structure. */
void c_geometry::reflect( int ix, int iy, int iz, int itx, int nop ) {
  int iti, i, nx, itagi, k;
  nec_float e1, e2, fnop, sam, cs, ss, xk, yk;
  
  np= n_segments;
  mp= m;
  m_ipsym=0;
  iti= itx;
  
  if ( ix >= 0) {
    if ( nop == 0)
    return;
  
    m_ipsym=1;
  
    /* reflect along z axis */
    if ( iz != 0)
    {
    m_ipsym=2;
  
    if ( n_segments > 0 )
    {
    /* Reallocate tags buffer */
    segment_tags.resize(2*n_segments + m);
    // segment_tags.resize((2*n_segments+m));
  
    /* Reallocate wire buffers */
    int new_size = 2*n_segments;
    x.resize(new_size);
    y.resize(new_size);
    z.resize(new_size);
    x2.resize(new_size);
    y2.resize(new_size);
    z2.resize(new_size);
    segment_radius.resize(new_size);
  
    for( i = 0; i < n_segments; i++ )
    {
    nx= i+ n_segments;
    e1= z[i];
    e2= z2[i];
  
    if ( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
    {
      nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
      nex->append(i+1);
      nex->append("LIES IN PLANE OF SYMMETRY");
      throw nex;
    }
  
    x[nx]= x[i];
    y[nx]= y[i];
    z[nx] = - e1;
    x2[nx]= x2[i];
    y2[nx]= y2[i];
    z2[nx] = - e2;
    itagi= segment_tags[i];
  
    if ( itagi == 0)
      segment_tags[nx]=0;
    if ( itagi != 0)
      segment_tags[nx]= itagi+ iti;
  
    segment_radius[nx]= segment_radius[i];
  
    } /* for( i = 0; i < n_segments; i++ ) */
  
    n_segments= n_segments*2;
    iti= iti*2;
  
    } /* if ( n_segments > 0) */
  
    if ( m > 0 )
    {
    /* Reallocate patch buffers */
    int new_size = 2*m;
    px.resize(new_size);
    py.resize(new_size);
    pz.resize(new_size);
    t1x.resize(new_size);
    t1y.resize(new_size);
    t1z.resize(new_size);
    t2x.resize(new_size);
    t2y.resize(new_size);
    t2z.resize(new_size);
    pbi.resize(new_size);
    psalp.resize(new_size);
  
    for( i = 0; i < m; i++ )
    {
      nx = i+m;
      if ( fabs(pz[i]) <= 1.0e-10)
      {
        nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--PATCH ");
        nex->append(i+1);
        nex->append("LIES IN PLANE OF SYMMETRY");
        throw nex;
      }
  
      px[nx]= px[i];
      py[nx]= py[i];
      pz[nx] = - pz[i];
      t1x[nx]= t1x[i];
      t1y[nx]= t1y[i];
      t1z[nx] = - t1z[i];
      t2x[nx]= t2x[i];
      t2y[nx]= t2y[i];
      t2z[nx] = - t2z[i];
      psalp[nx] = - psalp[i];
      pbi[nx]= pbi[i];
    }
  
    m= m*2;
  
    } /* if ( m >= m2) */
  
    } /* if ( iz != 0) */
  
    /* reflect along y axis */
    if ( iy != 0)  {
      if ( n_segments > 0)  {
      /* Reallocate tags buffer */
      segment_tags.resize(2*n_segments + m);
      // segment_tags.resize((2*n_segments+m));/*????*/
    
      /* Reallocate wire buffers */
      int new_size = 2*n_segments;
      x.resize(new_size);
      y.resize(new_size);
      z.resize(new_size);
      x2.resize(new_size);
      y2.resize(new_size);
      z2.resize(new_size);
      segment_radius.resize(new_size);
    
      for( i = 0; i < n_segments; i++ )  {
        nx= i+ n_segments;
        e1= y[i];
        e2= y2[i];
      
        if ( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )  {
          nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
          nex->append(i+1);
          nex->append("LIES IN PLANE OF SYMMETRY");
          throw nex;
        }
      
        x[nx]= x[i];
        y[nx] = - e1;
        z[nx]= z[i];
        x2[nx]= x2[i];
        y2[nx] = - e2;
        z2[nx]= z2[i];
        itagi= segment_tags[i];
      
        if ( itagi == 0)
          segment_tags[nx]=0;
        if ( itagi != 0)
          segment_tags[nx]= itagi+ iti;
      
        segment_radius[nx]= segment_radius[i];
      } /* for( i = n2-1; i < n_segments; i++ ) */
    
      n_segments= n_segments*2;
      iti= iti*2;
    } /* if ( n_segments >= n2) */
  
    if ( m > 0 )  {
      /* Reallocate patch buffers */
      int new_size = 2*m;
      px.resize(new_size);
      py.resize(new_size);
      pz.resize(new_size);
      t1x.resize(new_size);
      t1y.resize(new_size);
      t1z.resize(new_size);
      t2x.resize(new_size);
      t2y.resize(new_size);
      t2z.resize(new_size);
      pbi.resize(new_size);
      psalp.resize(new_size);
    
      for( i = 0; i < m; i++ )  {
        nx= i+m;
        if ( fabs( py[i]) <= 1.0e-10)  {
          nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--PATCH ");
          nex->append(i+1);
          nex->append("LIES IN PLANE OF SYMMETRY");
          throw nex;
        }
      
        px[nx]= px[i];
        py[nx] = - py[i];
        pz[nx]= pz[i];
        t1x[nx]= t1x[i];
        t1y[nx] = - t1y[i];
        t1z[nx]= t1z[i];
        t2x[nx]= t2x[i];
        t2y[nx] = - t2y[i];
        t2z[nx]= t2z[i];
        psalp[nx] = - psalp[i];
        pbi[nx]= pbi[i];
      } /* for( i = m2; i <= m; i++ ) */
    
      m= m*2;
    } /* if ( m >= m2) */
  
    } /* if ( iy != 0) */
  
    /* reflect along x axis */
    if ( ix == 0 )
      return;
  
    if ( n_segments > 0 )  {
      /* Reallocate tags buffer */
      segment_tags.resize(2*n_segments + m);
      // segment_tags.resize((2*n_segments+m));/*????*/
    
      /* Reallocate wire buffers */
      int new_size = 2*n_segments;
      x.resize(new_size);
      y.resize(new_size);
      z.resize(new_size);
      x2.resize(new_size);
      y2.resize(new_size);
      z2.resize(new_size);
      segment_radius.resize(new_size);
    
      for( i = 0; i < n_segments; i++ )  {
        nx= i+ n_segments;
        e1= x[i];
        e2= x2[i];
      
        if ( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
        {
          nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
          nex->append(i+1);
          nex->append("LIES IN PLANE OF SYMMETRY");
          throw nex;
        }
      
        x[nx] = - e1;
        y[nx]= y[i];
        z[nx]= z[i];
        x2[nx] = - e2;
        y2[nx]= y2[i];
        z2[nx]= z2[i];
        itagi= segment_tags[i];
      
        if ( itagi == 0)
        segment_tags[nx]=0;
        if ( itagi != 0)
        segment_tags[nx]= itagi+ iti;
      
        segment_radius[nx]= segment_radius[i];
      }
    
      n_segments= n_segments*2;
  
    } /* if ( n_segments > 0) */
  
    if ( m == 0 )
      return;
  
    /* Reallocate patch buffers */
    int new_size = 2*m;
    px.resize(new_size);
    py.resize(new_size);
    pz.resize(new_size);
    t1x.resize(new_size);
    t1y.resize(new_size);
    t1z.resize(new_size);
    t2x.resize(new_size);
    t2y.resize(new_size);
    t2z.resize(new_size);
    pbi.resize(new_size);
    psalp.resize(new_size);
  
    for( i = 0; i < m; i++ )  {
      nx= i+m;
      if ( fabs( px[i]) <= 1.0e-10)  {
        nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--PATCH ");
        nex->append(i+1);
        nex->append("LIES IN PLANE OF SYMMETRY");
        throw nex;
      }
    
      px[nx] = - px[i];
      py[nx]= py[i];
      pz[nx]= pz[i];
      t1x[nx] = - t1x[i];
      t1y[nx]= t1y[i];
      t1z[nx]= t1z[i];
      t2x[nx] = - t2x[i];
      t2y[nx]= t2y[i];
      t2z[nx]= t2z[i];
      psalp[nx] = - psalp[i];
      pbi[nx]= pbi[i];
    }
  
    m= m*2;
    return;
  } /* if ( ix >= 0) */
  
  /* reproduce structure with rotation to form cylindrical structure */
  fnop= (nec_float)nop;
  m_ipsym = -1;
  sam=two_pi() / fnop;
  cs= cos( sam);
  ss= sin( sam);
  
  if ( n_segments > 0)  {
    n_segments *= nop;
    nx= np;
  
    /* Reallocate tags buffer */
    segment_tags.resize(n_segments + m);
    //segment_tags.resize((n_segments+m));/*????*/
  
    /* Reallocate wire buffers */
    x.resize(n_segments);
    y.resize(n_segments);
    z.resize(n_segments);
    x2.resize(n_segments);
    y2.resize(n_segments);
    z2.resize(n_segments);
    segment_radius.resize(n_segments);
  
    for( i = nx; i < n_segments; i++ )  {
      k= i- np;
      xk= x[k];
      yk= y[k];
      x[i]= xk* cs- yk* ss;
      y[i]= xk* ss+ yk* cs;
      z[i]= z[k];
      xk= x2[k];
      yk= y2[k];
      x2[i]= xk* cs- yk* ss;
      y2[i]= xk* ss+ yk* cs;
      z2[i]= z2[k];
      segment_radius[i]= segment_radius[k];
      itagi= segment_tags[k];
    
      if ( itagi == 0)
        segment_tags[i]=0;
      if ( itagi != 0)
        segment_tags[i]= itagi+ iti;
    }
  
  } /* if ( n_segments >= n2) */
  
  if ( m == 0 )
    return;
  
  m *= nop;
  nx= mp;
  
  /* Reallocate patch buffers */
  px.resize(m);
  py.resize(m);
  pz.resize(m);
  t1x.resize(m);
  t1y.resize(m);
  t1z.resize(m);
  t2x.resize(m);
  t2y.resize(m);
  t2z.resize(m);
  pbi.resize(m);
  psalp.resize(m);
  
  for( i = nx; i < m; i++ )  {
    k = i-mp;
    xk= px[k];
    yk= py[k];
    px[i]= xk* cs- yk* ss;
    py[i]= xk* ss+ yk* cs;
    pz[i]= pz[k];
    xk= t1x[k];
    yk= t1y[k];
    t1x[i]= xk* cs- yk* ss;
    t1y[i]= xk* ss+ yk* cs;
    t1z[i]= t1z[k];
    xk= t2x[k];
    yk= t2y[k];
    t2x[i]= xk* cs- yk* ss;
    t2y[i]= xk* ss+ yk* cs;
    t2z[i]= t2z[k];
    psalp[i]= psalp[k];
    pbi[i]= pbi[k];
  
  } /* for( i = nx; i < m; i++ ) */
}
  
/*-----------------------------------------------------------------------*/

/*! \brief Scale all dimensions of a structure by a constant.*/
void c_geometry::scale( nec_float xw1 ) {
  // scale wires
  for (int i = 0; i < n_segments; i++) {
    x[i] = x[i]* xw1;
    y[i] = y[i]* xw1;
    z[i] = z[i]* xw1;
    x2[i] = x2[i]* xw1;
    y2[i] = y2[i]* xw1;
    z2[i] = z2[i]* xw1;
    segment_radius[i]= segment_radius[i] * xw1;
  }
    
  if ( m > 0) {
    // scale patches
    nec_float yw1= xw1* xw1;
    for (int i = 0; i < m; i++)
    {
      px[i]= px[i]* xw1;
      py[i]= py[i]* xw1;
      pz[i]= pz[i]* xw1;
      pbi[i]= pbi[i]* yw1;
    }
  } /* if ( m > 0) */
}

/*-----------------------------------------------------------------------*/

/* connect sets up segment connection data in arrays icon1 and */
/* icon2 by searching for segment ends that are in contact. */
void c_geometry::connect_segments( int ignd )
{
  nscon= -1;
  maxcon = 1;
  
  if (n_segments <= 1) {
    throw new nec_exception("GEOMETRY HAS ONE OR FEWER SEGMENTS. Please send bug report. This causes an error that we're trying to fix.");
  }
  
  if ( ignd != 0) {
    m_output->nec_printf( "\n\n     GROUND PLANE SPECIFIED." );
  
    if ( ignd > 0)
      m_output->nec_printf(
        "\n     WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL"
        " BE INTERPOLATED TO IMAGE IN GROUND PLANE.\n" );
  
    if ( m_ipsym == 2) {
      np=2* np;
      mp=2* mp;
    }
  
    if ( abs( m_ipsym) > 2 ) {
      np= n_segments;
      mp= m;
    }
  
    if ( np > n_segments) {
      throw new nec_exception("ERROR: NP > N IN c_geometry::connect_segments()" );
    }
  
    if ( (np == n_segments) && (mp == m) )
      m_ipsym=0;
  } /* if ( ignd != 0) */
  
  if ( n_segments != 0) {
    /* Allocate memory to connections */
    icon1.resize((n_segments+m));
    icon2.resize((n_segments+m));
  
    for (int i = 0; i < n_segments; i++ ) {
      int iz = i+1;
    
      nec_float zi1 = z[i];
      nec_float zi2 = z2[i];
      
      nec_3vector v1(x[i], y[i], z[i]);
      nec_3vector v2(x2[i], y2[i], z2[i]);
      nec_float slen = norm(v2 - v1) * SMIN;    
      
      /* determine connection data for end 1 of segment. */
      bool segment_on_ground = false;
      if ( ignd > 0) {
        if ( zi1 <= -slen) {
          nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
          nex->append(iz);
          nex->append("EXTENDS BELOW GROUND");
          throw nex;
        }
      
        if ( zi1 <= slen) {
          icon1[i]= iz;
          z[i]=0.;
          segment_on_ground = true;  
        } /* if ( zi1 <= slen) */
      } /* if ( ignd > 0) */
    
      if ( false == segment_on_ground ) {
        int ic= i;
        nec_float sep=0.0;
        for (int j = 1; j < n_segments; j++) {
          ic++;
          if ( ic >= n_segments)
            ic=0;
        
          nec_3vector vic(x[ic], y[ic], z[ic]);
          sep = normL1(v1 - vic);
          if ( sep <= slen) {
            icon1[i]= -(ic+1);
            break;
          }
        
          nec_3vector v2ic(x2[ic], y2[ic], z2[ic]);
          sep = normL1(v1 - v2ic);
          if ( sep <= slen) {
            icon1[i]= (ic+1);
            break;
          }
        
        } /* for( j = 1; j < n_segments; j++) */
      
        if ( ((iz > 0) || (icon1[i] <= PCHCON)) && (sep > slen) )
          icon1[i]=0;
      
      } /* if ( ! jump ) */
    
      /* determine connection data for end 2 of segment. */
      if ( (ignd > 0) || segment_on_ground ) {
        if ( zi2 <= -slen) {
          nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
          nex->append(iz);
          nex->append("EXTENDS BELOW GROUND");
          throw nex;
        }
      
        if ( zi2 <= slen) {
          if ( icon1[i] == iz ) {
            nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
            nex->append(iz);
            nex->append("LIES IN GROUND PLANE");
            throw nex;
          }
        
          icon2[i]= iz;
          z2[i]=0.;
          continue;
        } /* if ( zi2 <= slen) */  
      } /* if ( ignd > 0) */
    
      // re-initialize these vectors!
      v1 = nec_3vector(x[i], y[i], z[i]);
      v2 = nec_3vector(x2[i], y2[i], z2[i]);
      int ic= i;
      nec_float sep=0.0;
      for (int j = 1; j < n_segments; j++ ) {
        ic++;
        if ( ic >= n_segments)
          ic=0;
      
        nec_3vector vic(x[ic], y[ic], z[ic]);
        sep = normL1(v2 - vic);
        if (sep <= slen) {
          icon2[i]= (ic+1);
          break;
        }
      
        nec_3vector v2ic(x2[ic], y2[ic], z2[ic]);
        sep = normL1(v2 - v2ic);
        if (sep <= slen) {
          icon2[i]= -(ic+1);
          break;
        }
      } /* for( j = 1; j < n_segments; j++ ) */
    
      if ( ((iz > 0) || (icon2[i] <= PCHCON)) && (sep > slen) )
        icon2[i]=0;
    
    } /* for( i = 0; i < n_segments; i++ ) */
    

    /* find wire-surface connections for new patches */
    for (int ix=0; ix <m; ix++) {
//      DEBUG_TRACE("i: " << ix+1 << " ix: " << ix << " m: " << m);    
      nec_3vector vs(px[ix], py[ix], pz[ix]);
    
      for (int iseg = 0; iseg < n_segments; iseg++ ) {
        nec_3vector v1(x[iseg], y[iseg], z[iseg]);
        nec_3vector v2(x2[iseg], y2[iseg], z2[iseg]);
      
        /* for first end of segment */
        nec_float slen = normL1(v2 - v1) * SMIN;
        
        nec_float sep = normL1(v1 - vs);
        /* connection - divide patch into 4 patches at present array loc. */
        if ( sep <= slen) {
          icon1[iseg]=PCHCON + ix + 1;
          divide_patch(ix + 1);
          break;
        }
      
        sep = normL1(v2 - vs);
        
        if ( sep <= slen) {
          icon2[iseg]=PCHCON+ ix + 1;
          divide_patch(ix + 1);
          break;
        }
      }
    }  
  } /* if ( n_segments != 0) */
  
  m_output->nec_printf( "\n\n"
    "     TOTAL SEGMENTS USED: %d   SEGMENTS IN A"
    " SYMMETRIC CELL: %d   SYMMETRY FLAG: %d",
    n_segments, np, m_ipsym );
  
  if ( m > 0)
    m_output->nec_printf(  "\n"
    "       TOTAL PATCHES USED: %d   PATCHES"
    " IN A SYMMETRIC CELL: %d",  m, mp );
  
  
  if (0 == np + mp)
    throw new nec_exception("connect_segments Geometry has zero wires and zero patches.");
  
  int symmetry = (n_segments+m)/(np+mp); /* was iseg */
  if ( symmetry != 1)  {
    /*** may be error condition?? ***/
    if ( m_ipsym == 0 )  {
      nec_error_mode nem(*m_output);
      m_output->endl();
      m_output->line("ERROR: IPSYM=0 IN connect_segments()" );
      throw new nec_exception("ERROR: IPSYM=0 IN connect_segments()");
    }
  
    if ( m_ipsym < 0 )
      m_output->nec_printf(
        "\n  STRUCTURE HAS %d FOLD ROTATIONAL SYMMETRY\n", symmetry );
    else {
      int sym_planes = symmetry/2;
      if ( symmetry == 8)
        sym_planes=3;
      m_output->nec_printf(
        "\n  STRUCTURE HAS %d PLANES OF SYMMETRY\n", sym_planes );
    } /* if ( m_ipsym < 0 ) */
  
  } /* if ( symmetry == 1) */
    
  if ( n_segments == 0)
    return;
  
  /* Allocate to connection buffers */
  jco.resize(maxcon);
  
  int junction_counter = 0; // used just to print the junction number out if there are 3 or more segments
  bool header_printed = false; // Have we printed the header
  for (int j = 0; j < n_segments; j++ ) {
    int jx = j+1;
    int iend = -1;
    int jend = -1;
    int ix= icon1[j];
    int ic=1;
    jco[0]= -jx;
    nec_float xa = x[j];
    nec_float ya = y[j];
    nec_float za = z[j];
  
    while ( true )  {
      if ( (ix != 0) && (ix != (j+1)) && (ix <= PCHCON) )  {
        bool jump = false;
        // int nsflg = 0;  // will be set to 1 if the junction includes any new segments when NGF is in use.
        // NOTE nsflg is not used correctly as we don't use Numerical Greens Functions
        do {
          
          if ( ix == 0 ) {
            nec_exception* nex = new nec_exception("CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT: ");
            nex->append(ix);
            throw nex;
          }
        
          if ( ix < 0 )
            ix= -ix;
          else
            jend= -jend;
        
          jump = false;
        
          if ( ix == jx )
            break;
        
          if ( ix < jx ) {
            jump = true;
            break;
          }
        
          /* Record max. no. of connections */
          ic++;
          if ( ic >= maxcon ) {
            maxcon = ic+1;
            jco.resize(maxcon);
          }
          jco[ic-1]= ix* jend;
        
          // if ( ix > 0)
          //  nsflg=1;
        
          int ixx = ix-1;
          if ( jend != 1) {
            xa= xa+ x[ixx]; // dies here if n_segments == 1. ix is totally fried.
            ya= ya+ y[ixx];
            za= za+ z[ixx];
            ix= icon1[ixx];
          } else {
            xa= xa+ x2[ixx];
            ya= ya+ y2[ixx];
            za= za+ z2[ixx];
            ix= icon2[ixx];
          }
        }
        while( ix != 0 );
      
        if ( jump && (iend == 1) )
          break;
        else
          if ( jump ) {
            iend=1;
            jend=1;
            ix= icon2[j];
            ic=1;
            jco[0]= jx;
            xa= x2[j];
            ya= y2[j];
            za= z2[j];
            continue;
          }
      
        nec_float sep = (nec_float)ic;
        xa = xa / sep;
        ya = ya / sep;
        za = za / sep;
      
        for (int i = 0; i < ic; i++ ) {
          ix= jco[i];
          if ( ix <= 0) {
            ix = - ix;
            int ixx = ix-1; // TODO if ix == 0 we have a problem
            x[ixx]= xa;
            y[ixx]= ya;
            z[ixx]= za;
            continue;
          }
        
          int ixx = ix-1;
          x2[ixx]= xa;
          y2[ixx]= ya;
          z2[ixx]= za;  
        } /* for( i = 0; i < ic; i++ ) */
      
        if ( ic >= 3) {
          if ( false == header_printed ) {
            m_output->nec_printf( "\n\n"
              "    ---------- MULTIPLE WIRE JUNCTIONS ----------\n"
              "    JUNCTION  SEGMENTS (- FOR END 1, + FOR END 2)" );
            header_printed = true;
          }
        
          junction_counter++;
          m_output->nec_printf( "\n   %5d      ", junction_counter );
        
          for (int i = 1; i <= ic; i++ ) {
            m_output->nec_printf( "%5d", jco[i-1] );
            if ( !(i % 20) )
              m_output->nec_printf( "\n              " );
          }
        } /* if ( ic >= 3) */
      
        } /*if ( (ix != 0) && (ix != j) && (ix <= PCHCON) ) */
      
      if ( iend == 1)
        break;
    
      iend=1;
      jend=1;
      ix= icon2[j];
      ic=1;
      jco[0]= jx;
      xa= x2[j];
      ya= y2[j];
      za= z2[j];
    } /* while( true ) */
  } /* for( j = 0; j < n_segments; j++ ) */
  
  ax.resize(maxcon);
  bx.resize(maxcon);
  cx.resize(maxcon);
}

/* arc generates segment geometry data for an arc of segment_count segments */
void c_geometry::arc( int tag_id, int segment_count, nec_float rada,
    nec_float ang1, nec_float ang2, nec_float rad )
{
  int istart = n_segments;
  n_segments += segment_count;
  np= n_segments;
  mp= m;
  m_ipsym=0;
  
  if ( segment_count < 1)
    return;
  
  if ( fabs(ang2 - ang1) > 360.0)
  {
    throw new nec_exception("ERROR -- ARC ANGLE EXCEEDS 360 DEGREES");
  }
  
  /* Reallocate tags buffer */
  segment_tags.resize(n_segments+m);

  /* Reallocate wire buffers */
  x.resize(n_segments);
  y.resize(n_segments);
  z.resize(n_segments);
  x2.resize(n_segments);
  y2.resize(n_segments);
  z2.resize(n_segments);
  segment_radius.resize(n_segments);

  nec_float ang = degrees_to_rad(ang1);
  nec_float dang = degrees_to_rad(ang2- ang1) / segment_count;
  nec_float xs1= rada * cos(ang);
  nec_float zs1= rada * sin(ang);

  for(int i = istart; i < n_segments; i++ )
  {
    ang += dang;
    nec_float xs2 = rada * cos(ang);
    nec_float zs2 = rada * sin(ang);
    x[i]= xs1;
    y[i]=0.;
    z[i]= zs1;
    x2[i]= xs2;
    y2[i]=0.;
    z2[i]= zs2;
    xs1= xs2;
    zs1= zs2;
    segment_radius[i]= rad;
    segment_tags[i]= tag_id;
  } /* for( i = ist; i < n_segments; i++ ) */
}


/*-----------------------------------------------------------------------*/

void c_geometry::sp_card(int ns,
      nec_float in_x1, nec_float in_y1, nec_float in_z1,
      nec_float in_x2, nec_float in_y2, nec_float in_z2)
{
  const char ipt[4] = { 'P', 'R', 'T', 'Q' };
  this->patch_type = ns;
  m_output->nec_printf( "\n"
          " %5d%c %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
          m+1, ipt[ns], in_x1, in_y1, in_z1, in_x2, in_y2, in_z2 );
  switch (ns) {
    case 0: { // Arbitrary Shape
      nec_float elevation = degrees_to_rad(in_x2);
      nec_float azimuth = degrees_to_rad(in_y2);
      patch( 0, 0, in_x1, in_y1, in_z1, elevation, azimuth, in_z2, 0, 0, 0, 0, 0, 0);
      break;
    }
      
    case 1: // Rectangular (require SC card)
    case 2: // triangular, (require SC card)
    case 3: // quadrilateral (require SC card)
      this->patch_x1 = nec_3vector(in_x1, in_y1, in_z1);
      this->patch_x2 = nec_3vector(in_x2, in_y2, in_z2);
      this->patch_type = ns;
      break;
      
    default:
      throw new nec_exception("PATCH DATA ERROR ns ",ns );
  }
  _prev_sc = false;
}


void c_geometry::gx_card(int card_int_1, int card_int_2) {       
  const char ifx[2] = {'*', 'X'}, ify[2]={'*','Y'}, ifz[2]={'*','Z'};
  
  int iy= card_int_2/10;
  int iz= card_int_2- iy*10;
  int ix= iy/10;
  iy= iy- ix*10;

  if ( ix != 0)   ix=1;
  if ( iy != 0)   iy=1;
  if ( iz != 0)   iz=1;

  m_output->nec_printf(
          "\n  STRUCTURE REFLECTED ALONG THE AXES %c %c %c"
          " - TAGS INCREMENTED BY %d\n",
          ifx[ix], ify[iy], ifz[iz], card_int_1 );

  reflect( ix, iy, iz, card_int_1, card_int_2);
}
/*
 * For the rectangular or quadrilateral options, multiple SC cards may follow a SP card 
 * to specify a string of patches. The parameters on the second or subsequent SC card 
 * specify corner 3 for a rectangle or corners 3 and 4 for a quadrilateral, while 
 * corners 3 and 4 of the previous patch become corners 2 and 1, respectively, of 
 * the new patch. The integer I2 on the second or subsequent SC card specifies 
 * the new patch shape and must be 1 for rectangular shape or 3 for quadrilateral 
 * shape. On the first SC card after SP, I2 has no effect. 
 * 
 * Rectangular or quadrilateral patches may be intermixed, but triangular or arbitrary
 * shapes are not allowed in a string of linked patches. 
 */
void c_geometry::sc_card(int i2,
      nec_float x3, nec_float y3, nec_float z3,
      nec_float x4, nec_float y4, nec_float z4)
{
  if (_prev_sc) {
    sc_multiple_card(i2, x3, y3, z3, x4, y4, z4);
    return;
  }
    
  DEBUG_TRACE("sc_card(" << i2 << ")");
  
  m_output->nec_printf( "\n"
        "      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
        x3, y3, z3, x4, y4, z4 );
  
  this->patch_x3 = nec_3vector(x3, y3, z3);
  this->patch_x4 = nec_3vector(x4, y4, z4);

  switch (this->patch_type) {
    case 0:  // Arbritrary
      throw new nec_exception("PATCH DATA ERROR: SC CARD FOR ARBITRARY PATCH TYPE (NS=0) ");
      break;
    
    case 1:  // Rectangular
      this->patch_x4 = this->patch_x1 + this->patch_x3 - this->patch_x2;
      break;
      
    case 2:  // Triangular
      break;
      
    case 3:  // Quadrilateral
      break;
      
    default:
      throw new nec_exception("PATCH DATA ERROR ns ", this->patch_type );
  }
  patch( this->patch_type, i2, 
          this->patch_x1(0), this->patch_x1(1), this->patch_x1(2), 
          this->patch_x2(0), this->patch_x2(1), this->patch_x2(2),
          this->patch_x3(0), this->patch_x3(1), this->patch_x3(2),
          this->patch_x4(0), this->patch_x4(1), this->patch_x4(2));
  
  _prev_sc = true;
}

/*
 * For the rectangular or quadrilateral options, multiple SC cards may follow a SP card 
 * to specify a string of patches. The parameters on the second or subsequent SC card 
 * specify corner 3 for a rectangle or corners 3 and 4 for a quadrilateral, while 
 * corners 3 and 4 of the previous patch become corners 2 and 1, respectively, of 
 * the new patch. The integer I2 on the second or subsequent SC card specifies 
 * the new patch shape and must be 1 for rectangular shape or 3 for quadrilateral 
 * shape. On the first SC card after SP, I2 has no effect. 
 * 
 * Rectangular or quadrilateral patches may be intermixed, but triangular or arbitrary
 * shapes are not allowed in a string of linked patches. 
 */
void c_geometry::sc_multiple_card(int i2,
      nec_float x3, nec_float y3, nec_float z3,
      nec_float x4, nec_float y4, nec_float z4)
{
  DEBUG_TRACE("sc_multiple_card(" << i2 << ")");
  
  const char ipt[4] = { 'P', 'R', 'T', 'Q' };
  switch (i2) {
    case 0:  // Arbritrary
      throw new nec_exception("PATCH DATA ERROR: MULTIPLE SC CARDS FOR ARBITRARY PATCH TYPE (NS=0) ");
      break;
    case 2:  // Triangular
      throw new nec_exception("PATCH DATA ERROR: MULTIPLE SC CARDS FOR TRIANGULAR PATCH (NS=2) ");
      break;
    
    case 1:  // Rectangular
      this->patch_x1 = this->patch_x4;
      this->patch_x2 = this->patch_x3;
      this->patch_x3 = nec_3vector(x3, y3, z3);
      this->patch_x4 = this->patch_x1 + this->patch_x3 - this->patch_x2;
      break;
      
    case 3:  // Quadrilateral
      this->patch_x1 = this->patch_x4;
      this->patch_x2 = this->patch_x3;
      this->patch_x3 = nec_3vector(x3, y3, z3);
      this->patch_x4 = nec_3vector(x4, y4, z4);
      break;
      
    default:
      throw new nec_exception("PATCH DATA ERROR i2 = ", i2 );
  }
  
  patch( this->patch_type, i2, 
          this->patch_x1(0), this->patch_x1(1), this->patch_x1(2), 
          this->patch_x2(0), this->patch_x2(1), this->patch_x2(2),
          this->patch_x3(0), this->patch_x3(1), this->patch_x3(2),
          this->patch_x4(0), this->patch_x4(1), this->patch_x4(2));
  
  m_output->nec_printf( "\n"
          " %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
          this->patch_type, ipt[i2], 
          this->patch_x1(0), this->patch_x1(1), this->patch_x1(2),
          this->patch_x2(0), this->patch_x2(1), this->patch_x2(2));

  m_output->nec_printf( "\n"
          "      %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f",
          x3, y3, z3, x4, y4, z4 );

  _prev_sc = true;
}


/*! \brief patch generates and modifies patch geometry data.
*/
void c_geometry::patch( int nx, int ny,
    nec_float ax1, nec_float ay1, nec_float az1,
    nec_float ax2, nec_float ay2, nec_float az2,
    nec_float ax3, nec_float ay3, nec_float az3,
    nec_float ax4, nec_float ay4, nec_float az4 )
{
  int mi, ntp, iy, ix;
  nec_float s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., xst=0.;
  nec_float znv, xnv, ynv, xa, xn2, yn2, zn2, salpn, xs, ys, zs, xt, yt, zt;
  
  /* new patches.  for nx=0, ny=1,2,3,4 patch is (respectively) */;
  /* arbitrary, rectagular, triangular, or quadrilateral. */
  /* for nx and ny  > 0 a rectangular surface is produced with */
  /* nx by ny rectangular patches. */
  
  m++;
  mi= m-1;
  
  /* Reallocate patch buffers */
  px.resize(m);
  py.resize(m);
  pz.resize(m);
  t1x.resize(m);
  t1y.resize(m);
  t1z.resize(m);
  t2x.resize(m);
  t2y.resize(m);
  t2z.resize(m);
  pbi.resize(m);
  psalp.resize(m);
  
  if ( nx > 0)
    ntp=2;
  else
    ntp= ny;
  
  if ( ntp <= 1)
  {
    px[mi]= ax1;
    py[mi]= ay1;
    pz[mi]= az1;
    pbi[mi]= az2;
    znv= cos( ax2);
    xnv= znv* cos( ay2);
    ynv= znv* sin( ay2);
    znv= sin( ax2);
    xa= sqrt( xnv* xnv+ ynv* ynv);
  
    if ( xa >= 1.0e-6)
    {
      t1x[mi] = - ynv/ xa;
      t1y[mi]= xnv/ xa;
      t1z[mi]=0.;
    }
    else
    {
      t1x[mi]=1.;
      t1y[mi]=0.;
      t1z[mi]=0.;
    }
  
  } /* if ( ntp <= 1) */
  else
  {
    s1x= ax2- ax1;
    s1y= ay2- ay1;
    s1z= az2- az1;
    s2x= ax3- ax2;
    s2y= ay3- ay2;
    s2z= az3- az2;
  
    if ( nx != 0)
    {
      s1x= s1x/ nx;
      s1y= s1y/ nx;
      s1z= s1z/ nx;
      s2x= s2x/ ny;
      s2y= s2y/ ny;
      s2z= s2z/ ny;
    }
  
    xnv= s1y* s2z- s1z* s2y;
    ynv= s1z* s2x- s1x* s2z;
    znv= s1x* s2y- s1y* s2x;
    xa= sqrt( xnv* xnv+ ynv* ynv+ znv* znv);
    xnv= xnv/ xa;
    ynv= ynv/ xa;
    znv= znv/ xa;
    xst= sqrt( s1x* s1x+ s1y* s1y+ s1z* s1z);
    t1x[mi]= s1x/ xst;
    t1y[mi]= s1y/ xst;
    t1z[mi]= s1z/ xst;
  
    if ( ntp <= 2)
    {
      px[mi]= ax1+.5*( s1x+ s2x);
      py[mi]= ay1+.5*( s1y+ s2y);
      pz[mi]= az1+.5*( s1z+ s2z);
      pbi[mi]= xa;
    }
    else
    {
      if ( ntp != 4)
      {
        px[mi]=( ax1+ ax2+ ax3)/3.;
        py[mi]=( ay1+ ay2+ ay3)/3.;
        pz[mi]=( az1+ az2+ az3)/3.;
        pbi[mi]=.5* xa;
      }
      else
      {
        s1x= ax3- ax1;
        s1y= ay3- ay1;
        s1z= az3- az1;
        s2x= ax4- ax1;
        s2y= ay4- ay1;
        s2z= az4- az1;
        xn2= s1y* s2z- s1z* s2y;
        yn2= s1z* s2x- s1x* s2z;
        zn2= s1x* s2y- s1y* s2x;
        xst= sqrt( xn2* xn2+ yn2* yn2+ zn2* zn2);
        salpn=1./(3.*( xa+ xst));
        px[mi]=( xa*( ax1+ ax2+ ax3)+ xst*( ax1+ ax3+ ax4))* salpn;
        py[mi]=( xa*( ay1+ ay2+ ay3)+ xst*( ay1+ ay3+ ay4))* salpn;
        pz[mi]=( xa*( az1+ az2+ az3)+ xst*( az1+ az3+ az4))* salpn;
        pbi[mi]=.5*( xa+ xst);
        s1x=( xnv* xn2+ ynv* yn2+ znv* zn2)/ xst;
      
        if ( s1x <= 0.9998)
        {
          throw new nec_exception("ERROR -- CORNERS OF QUADRILATERAL PATCH DO NOT LIE IN A PLANE");
        }
      } /* if ( ntp != 4) */
    } /* if ( ntp <= 2) */
  
  } /* if ( ntp <= 1) */
  
  t2x[mi]= ynv* t1z[mi]- znv* t1y[mi];
  t2y[mi]= znv* t1x[mi]- xnv* t1z[mi];
  t2z[mi]= xnv* t1y[mi]- ynv* t1x[mi];
  psalp[mi]=1.;
  
  if ( nx != 0)
  {
    m += nx*ny-1;
  
    /* Reallocate patch buffers */
    px.resize(m);
    py.resize(m);
    pz.resize(m);
    t1x.resize(m);
    t1y.resize(m);
    t1z.resize(m);
    t2x.resize(m);
    t2y.resize(m);
    t2z.resize(m);
    pbi.resize(m);
    psalp.resize(m);
  
    xn2= px[mi]- s1x- s2x;
    yn2= py[mi]- s1y- s2y;
    zn2= pz[mi]- s1z- s2z;
    xs= t1x[mi];
    ys= t1y[mi];
    zs= t1z[mi];
    xt= t2x[mi];
    yt= t2y[mi];
    zt= t2z[mi];
  
    for( iy = 0; iy < ny; iy++ )
    {
      xn2 += s2x;
      yn2 += s2y;
      zn2 += s2z;
    
      for( ix = 1; ix <= nx; ix++ )
      {
        xst= (nec_float)ix;
        px[mi]= xn2+ xst* s1x;
        py[mi]= yn2+ xst* s1y;
        pz[mi]= zn2+ xst* s1z;
        pbi[mi]= xa;
        psalp[mi]=1.;
        t1x[mi]= xs;
        t1y[mi]= ys;
        t1z[mi]= zs;
        t2x[mi]= xt;
        t2y[mi]= yt;
        t2z[mi]= zt;
        mi++;
      } /* for( ix = 0; ix < nx; ix++ ) */
    } /* for( iy = 0; iy < ny; iy++ ) */
  
  } /* if ( nx != 0) */
  
  m_ipsym=0;
  np= n_segments;
  mp= m;
}


/*!\brief Divide a patch into four (was subph).
Used when a patch is connected to a wire. The patch
nx is  divided into 4 patches that become nx, nx+1, nx+2, nx+3. The other
patches are shifted to make room for the three new patches.

\param nx The index of the patch to divide (starting at 1)
*/
void c_geometry::divide_patch(int nx)

{  
  m += 3;

  px.resize(m);
  py.resize(m);
  pz.resize(m);
  t1x.resize(m);
  t1y.resize(m);
  t1z.resize(m);
  t2x.resize(m);
  t2y.resize(m);
  t2z.resize(m);
  pbi.resize(m);
  psalp.resize(m);
  
  /* Shift patches to make room for new ones */
  for (int iy = m-1; iy > nx; iy--)
  {
    int old_index = iy - 3;
    px[iy] = px[old_index];
    py[iy] = py[old_index];
    pz[iy] = pz[old_index];
    pbi[iy] = pbi[old_index];
    psalp[iy] = psalp[old_index];
    t1x[iy] = t1x[old_index];
    t1y[iy] = t1y[old_index];
    t1z[iy] = t1z[old_index];
    t2x[iy] = t2x[old_index];
    t2y[iy] = t2y[old_index];
    t2z[iy] = t2z[old_index];
  }
  
  /* divide patch for connection */
  int patch_index = nx-1;
  nec_float xs = px[patch_index];
  nec_float ys = py[patch_index];
  nec_float zs = pz[patch_index];
  nec_float xa = pbi[patch_index]/4.;
  nec_float xst = sqrt(xa)/2.;
  nec_float s1x = t1x[patch_index];
  nec_float s1y = t1y[patch_index];
  nec_float s1z = t1z[patch_index];
  nec_float s2x = t2x[patch_index];
  nec_float s2y = t2y[patch_index];
  nec_float s2z = t2z[patch_index];
  nec_float saln = psalp[patch_index];
  nec_float xt = xst;
  nec_float yt = xst;
  
  int new_index = patch_index;
  
  /* Generate the four new patches */
  for (int ix = 1; ix <= 4; ix++ )
  {
    px[new_index]= xs+ xt* s1x+ yt* s2x;
    py[new_index]= ys+ xt* s1y+ yt* s2y;
    pz[new_index]= zs+ xt* s1z+ yt* s2z;
    pbi[new_index]= xa;
    t1x[new_index]= s1x;
    t1y[new_index]= s1y;
    t1z[new_index]= s1z;
    t2x[new_index]= s2x;
    t2y[new_index]= s2y;
    t2z[new_index]= s2z;
    psalp[new_index]= saln;
  
    if (2 == ix)
      yt = -yt;
  
    if ( (ix == 1) || (ix == 3) )
      xt = -xt;
  
    new_index++;
  }
  
  /* Readjust the mp patch index to account for the added patches */
  if (nx <= mp)
    mp += 3;
}

/*-----------------------------------------------------------------------

  Read Geometry Data from a Card

-------------------------------------------------------------------------*/

void c_geometry::read_geometry_card(FILE* input_fp,  char *gm,
  int *in_i1, int *in_i2,
  nec_float *in_x1, nec_float *in_y1, nec_float *in_z1,
  nec_float *in_x2, nec_float *in_y2, nec_float *in_z2,
  nec_float *in_rad )
{
  char line_buf[134];
  int i, line_idx;
  int n_integer_params = 2, n_float_params = 7;
  int integer_params[2] = { 0, 0 };
  nec_float real_params[7] = { 0., 0., 0., 0., 0., 0., 0. };
  
  /* read a line from input file */
  load_line( line_buf, input_fp );
  
  /* get line length */
  int line_length = (int)strlen( line_buf );
  
  /* abort if card's mnemonic too short or missing */
  if ( line_length < 2 )
  {
    nec_exception* nex = new nec_exception("GEOMETRY DATA CARD ERROR:");
    nex->append(" CARD'S MNEMONIC CODE TOO SHORT OR MISSING.");
    throw nex;
  }
  
  /* extract card's mnemonic code */
  strncpy( gm, line_buf, 2 );
  gm[2] = '\0';
  
  /* Exit if "XT" command read (for testing) */
  if ( strcmp( gm, "XT" ) == 0 )
  {
      nec_exception* nex = new nec_exception("Exiting after an \"XT\" command in read_geometry_card()");
      throw nex;
  }
  
  /* Return if only mnemonic on card */
  if ( line_length == 2 )
  {
    *in_i1 = *in_i2 = 0;
    *in_x1 = *in_y1 = *in_z1 = *in_x2 = *in_y2 = *in_z2 = *in_rad = 0.0;
    return;
  }
  
  /* read integers from line */
  line_idx = 1;
  for( i = 0; i < n_integer_params; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
      (line_buf[  line_idx] >  '9')) &&
      (line_buf[  line_idx] != '+')  &&
      (line_buf[  line_idx] != '-') )
    if ( line_buf[line_idx] == '\0' )
    {
      *in_i1= integer_params[0];
      *in_i2= integer_params[1];
      *in_x1= real_params[0];
      *in_y1= real_params[1];
      *in_z1= real_params[2];
      *in_x2= real_params[3];
      *in_y2= real_params[4];
      *in_z2= real_params[5];
      *in_rad= real_params[6];
      return;
    }
  
    /* read an integer from line */
    integer_params[i] = atoi( &line_buf[line_idx] );
  
    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while(   (line_buf[++line_idx] != ' ') &&
        (line_buf[  line_idx] != ',') &&
        (line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if (   ((line_buf[line_idx] <  '0')  ||
          (line_buf[line_idx] >  '9')) &&
          (line_buf[line_idx] != '+')  &&
          (line_buf[line_idx] != '-') )
      {
        nec_stop(
          "GEOMETRY DATA CARD \"%s\" ERROR:"
          "\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
          gm, line_buf[line_idx], (line_idx+1)  );
      }
    } /* while( (line_buff[++line_idx] ... */
  
    /* Return on end of line */
    if ( line_buf[line_idx] == '\0' )
    {
      *in_i1= integer_params[0];
      *in_i2= integer_params[1];
      *in_x1= real_params[0];
      *in_y1= real_params[1];
      *in_z1= real_params[2];
      *in_x2= real_params[3];
      *in_y2= real_params[4];
      *in_z2= real_params[5];
      *in_rad= real_params[6];
      return;
    }
  
  } /* for( i = 0; i < n_integer_params; i++ ) */
  
  /* read nec_floats from line */
  for( i = 0; i < n_float_params; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
      (line_buf[  line_idx] >  '9')) &&
      (line_buf[  line_idx] != '+')  &&
      (line_buf[  line_idx] != '-')  &&
      (line_buf[  line_idx] != '.') )
    if ( line_buf[line_idx] == '\0' )
    {
      *in_i1= integer_params[0];
      *in_i2= integer_params[1];
      *in_x1= real_params[0];
      *in_y1= real_params[1];
      *in_z1= real_params[2];
      *in_x2= real_params[3];
      *in_y2= real_params[4];
      *in_z2= real_params[5];
      *in_rad= real_params[6];
      return;
    }
  
    /* read a nec_float from line */
    real_params[i] = atof( &line_buf[line_idx] );
  
    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while(  (line_buf[++line_idx] != ' ') &&
        (line_buf[  line_idx] != ',') &&
        (line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if ( ((line_buf[line_idx] <  '0')  ||
        (line_buf[line_idx] >  '9')) &&
        (line_buf[line_idx] != '.')  &&
        (line_buf[line_idx] != '+')  &&
        (line_buf[line_idx] != '-')  &&
        (line_buf[line_idx] != 'E')  &&
        (line_buf[line_idx] != 'e') )
      {
        nec_stop(
          "\n  GEOMETRY DATA CARD \"%s\" ERROR:"
          "\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d.\n",
          gm, line_buf[line_idx], (line_idx+1) );
      }
    } /* while( (line_buff[++line_idx] ... */
  
    /* Return on end of line */
    if ( line_buf[line_idx] == '\0' )
    {
      *in_i1= integer_params[0];
      *in_i2= integer_params[1];
      *in_x1= real_params[0];
      *in_y1= real_params[1];
      *in_z1= real_params[2];
      *in_x2= real_params[3];
      *in_y2= real_params[4];
      *in_z2= real_params[5];
      *in_rad= real_params[6];
      return;
    }
  
  } /* for( i = 0; i < n_float_params; i++ ) */
  
  *in_i1  = integer_params[0];
  *in_i2  = integer_params[1];
  *in_x1  = real_params[0];
  *in_y1  = real_params[1];
  *in_z1  = real_params[2];
  *in_x2  = real_params[3];
  *in_y2  = real_params[4];
  *in_z2  = real_params[5];
  *in_rad = real_params[6];
}



/* compute basis function i */
void c_geometry::tbf( int i, int icap )
{
  int jcoxx, njun1=0, njun2, jsnop, jsnox;
  nec_float sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;
  
  jsno=0;
  
  nec_float pp = 0.0;
  int ix = i-1;
  int jcox = icon1[ix];
  
  if ( jcox > PCHCON)
    jcox= i;
  
  int jend = -1;
  int iend = -1;
  
  nec_float _sig = -1.0;
  do
  {
    if ( jcox != 0 )
    {
      if ( jcox < 0 )
        jcox = - jcox;
      else
      {
        _sig = -_sig;
        jend = - jend;
      }
    
      jcoxx = jcox-1;
      jsno++;
      jsnox = jsno-1;
      jco[jsnox]= jcox;
      nec_float d = pi() * segment_length[jcoxx];
      sdh = sin(d);
      cdh = cos(d);
      sd = 2.0 * sdh * cdh;
    
      if ( d <= 0.015)
      {
        omc=4.* d* d;
        omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
      }
      else
        omc=1.- cdh* cdh+ sdh* sdh;
    
      aj=1./( log(1./( pi()* segment_radius[jcoxx]))-.577215664);
      pp= pp- omc/ sd* aj;
      ax[jsnox]= aj/ sd* _sig;
      bx[jsnox]= aj/(2.* cdh);
      cx[jsnox] = - aj/(2.* sdh)* _sig;
    
      if ( jcox != i)
      {
        if ( jend == 1)
          jcox= icon2[jcoxx];
        else
          jcox= icon1[jcoxx];
      
        if ( abs(jcox) != i )
        {
          if ( jcox != 0 )
            continue;
          else
          {
            nec_exception* nex = new nec_exception("TBF - SEGMENT CONNECTION ERROR FOR SEGMENT ");
            nex->append(i);
            throw nex;
          }
        }
      } /* if ( jcox != i) */
      else
      bx[jsnox] = - bx[jsnox];
    
      if ( iend == 1)
        break;
    } /* if ( jcox != 0 ) */
  
    pm = - pp;
    pp=0.;
    njun1= jsno;
  
    jcox= icon2[ix];
    if ( jcox > PCHCON)
      jcox= i;
  
    jend=1;
    iend=1;
    _sig = -1.;
  
  } /* do */
  while( jcox != 0 );
  
  njun2= jsno- njun1;
  jsnop= jsno;
  jco[jsnop]= i;
  
  nec_float d = pi()* segment_length[ix];
  sdh = sin(d);
  cdh = cos(d);
  sd = 2.0 * sdh * cdh;
  cd= cdh*cdh - sdh*sdh;
  
  if ( d <= 0.015)
  {
    omc = 4.0* d*d;
    omc = ((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
  }
  else
    omc = 1.0 - cd;
  
  ap=1./( log(1./( pi()* segment_radius[ix]))-.577215664);
  aj= ap;
  
  if ( njun1 == 0)
  {
    if ( njun2 == 0)
    {
      bx[jsnop]=0.;
    
      if ( icap == 0)
        xxi=0.;
      else
      {
        qp= pi()* segment_radius[ix];
        xxi= qp* qp;
        xxi= qp*(1.-.5* xxi)/(1.- xxi);
      }
    
      cx[jsnop]=1./( cdh- xxi* sdh);
      jsno= jsnop+1;
      ax[jsnop] = -1.;
      return;
    } /* if ( njun2 == 0) */
  
    if ( icap == 0)
      xxi=0.;
    else
    {
      qp= pi()* segment_radius[ix];
      xxi= qp* qp;
      xxi= qp*(1.-.5* xxi)/(1.- xxi);
    }
  
    qp = -( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));
    d= cd- xxi* sd;
    bx[jsnop]=( sdh+ ap* qp*( cdh- xxi* sdh))/ d;
    cx[jsnop]=( cdh+ ap* qp*( sdh+ xxi* cdh))/ d;
  
    for( iend = 0; iend < njun2; iend++ )
    {
      ax[iend] = -ax[iend]* qp;
      bx[iend]= bx[iend]* qp;
      cx[iend] = - cx[iend]* qp;
    }
  
    jsno= jsnop+1;
    ax[jsnop] = -1.;
    return;
  } /* if ( njun1 == 0) */
  
  if ( njun2 == 0)
  {
    if ( icap == 0)
      xxi=0.;
    else
    {
      qm= pi()* segment_radius[ix];
      xxi= qm* qm;
      xxi= qm*(1.-.5* xxi)/(1.- xxi);
    }
  
    qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));
    d= cd- xxi* sd;
    bx[jsnop]=( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
    cx[jsnop]=( cdh- aj* qm*( sdh+ xxi* cdh))/ d;
  
    for( iend = 0; iend < njun1; iend++ )
    {
      ax[iend]= ax[iend]* qm;
      bx[iend]= bx[iend]* qm;
      cx[iend]= cx[iend]* qm;
    }
  
    jsno= jsnop+1;
    ax[jsnop] = -1.;
    return;
  
  } /* if ( njun2 == 0) */
  
  qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
  qm=( ap* omc- pp* sd)/ qp;
  qp = -( aj* omc+ pm* sd)/ qp;
  bx[jsnop]=( aj* qm+ ap* qp)* sdh/ sd;
  cx[jsnop]=( aj* qm- ap* qp)* cdh/ sd;
  
  for( iend = 0; iend < njun1; iend++ )
  {
    ax[iend]= ax[iend]* qm;
    bx[iend]= bx[iend]* qm;
    cx[iend]= cx[iend]* qm;
  }
  
  jend= njun1;
  for( iend = jend; iend < jsno; iend++ )
  {
    ax[iend] = - ax[iend]* qp;
    bx[iend]= bx[iend]* qp;
    cx[iend] = - cx[iend]* qp;
  }
  
  jsno= jsnop+1;
  ax[jsnop] = -1.;
}



/* compute the components of all basis functions on segment j */
void c_geometry::trio( int j )  {
  int jcox, jcoxx, jsnox, jx, jend=0, iend=0;
  
  jsno=0;
  jx = j-1;
  jcox= icon1[jx];
  jcoxx = jcox-1;
  
  if ( jcox <= PCHCON)  {
    jend = -1;
    iend = -1;
  }
  
  if ( (jcox == 0) || (jcox > PCHCON) )  {
    jcox= icon2[jx];
    jcoxx = jcox-1;
  
    if ( jcox <= PCHCON)  {
      jend=1;
      iend=1;
    }
  
    if ( jcox == 0 || (jcox > PCHCON) )  {
      jsnox = jsno;
      jsno++;
    
      /* Allocate to connections buffers */
      if ( jsno >= maxcon )  {
        maxcon = jsno +1;
        jco.resize(maxcon);
        ax.resize(maxcon);
        bx.resize(maxcon);
        cx.resize(maxcon);
      }
    
      sbf( j, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
      jco[jsnox]= j;
      return;
    }
  } /* if ( (jcox == 0) || (jcox > PCHCON) ) */
  
  do  {
    if ( jcox < 0 )
      jcox = - jcox;
    else
      jend = - jend;
    jcoxx = jcox-1;
  
    if ( jcox != j)  {
      jsnox = jsno;
      jsno++;
    
      /* Allocate to connections buffers */
      if ( jsno >= maxcon )  {
        maxcon = jsno +1;
        jco.resize(maxcon );
        ax.resize(maxcon);
        bx.resize(maxcon);
        cx.resize(maxcon);
      }
    
      sbf( jcox, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
      jco[jsnox]= jcox;
    
      if ( jend != 1)
        jcox= icon1[jcoxx];
      else
        jcox= icon2[jcoxx];
    
      if ( jcox == 0 )  {
        nec_exception* nex = new nec_exception("TRIO - SEGMENT CONNENTION ERROR FOR SEGMENT ");
        nex->append(j);
        throw nex;
      }
      else
        continue;
    } /* if ( jcox != j) */
  
    if ( iend == 1)
      break;
  
    jcox= icon2[jx];
  
    if ( jcox > PCHCON)
      break;
  
    jend=1;
    iend=1;
  } /* do */
  while( jcox != 0 );
  
  jsnox = jsno;
  jsno++;
  
  /* Allocate to connections buffers */
  if ( jsno >= maxcon )  {
    maxcon = jsno +1;
    jco.resize(maxcon );
    ax.resize(maxcon);
    bx.resize(maxcon);
    cx.resize(maxcon);
  }
  
  sbf( j, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
  jco[jsnox]= j;
}



/*! \brief To evaluate the current expansion function associated with a given segment, returning only that portion on a particular segment.
\param i The segment on which the expansion function is centered.
\param is The segment for which the function coefficients A_j, B_j and C_j are requested
\param aa The return value for the function coefficient A_j
\param bb The return value for the function coefficient B_j
\param cc The return value for the function coefficient C_j

SBF is very similar to TBF. Both routines evaluate the current expansion functions. However, while TBF stores the coefficients for each segment on which a given expansion function is non-zero, SBF returns the coefficients ofr only a single specified segment.

Refer to TBF for a discussion of the coding and variables. One additional variable in SBF -- june -- is set to -1 or +1 if \param is is found connected to end1 or end2 respectively of segment i. If I == IS and segment i is not connected to a surface of ground plane, then june is set to 0.
*/
void c_geometry::sbf( int i, int is, nec_float *aa, nec_float *bb, nec_float *cc )
{
  int local_jsno; // this parameter is renamed because it shadows the member variable of the same name
  int ix, june, jcox, jcoxx, jend, iend, njun1=0, njun2;
  nec_float d, sig, pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;
  
  *aa=0.;
  *bb=0.;
  *cc=0.;
  june=0;
  local_jsno=0;
  pp=0.;
  ix=i-1;
  
  jcox= icon1[ix];
  if ( jcox > PCHCON)
    jcox= i;
  jcoxx = jcox-1;
  
  jend = -1;
  iend = -1;
  sig = -1.;
  
  do  {
    // DEBUG_TRACE("c_geometry::sbf(" << i << "," << is << "): " << jcox);
    if ( jcox != 0 )  {
      if ( jcox < 0 )
        jcox = - jcox;
      else  {
        sig = - sig;
        jend = - jend;
      }
    
      jcoxx = jcox-1;
      local_jsno++;
      d= pi()* segment_length[jcoxx];
      sdh= sin( d);
      cdh= cos( d);
      sd=2.* sdh* cdh;
    
      if ( d <= 0.015)  {
        omc=4.* d* d;
        omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
      }
      else
        omc=1.- cdh* cdh+ sdh* sdh;
    
      aj=1./( log(1./( pi()* segment_radius[jcoxx]))-.577215664);
      pp -= omc/ sd* aj;
    
      if ( jcox == is)  {
        *aa= aj/ sd* sig;
        *bb= aj/(2.* cdh);
        *cc = - aj/(2.* sdh)* sig;
        june= iend;
      }
    
      if ( jcox != i )  {
        if ( jend != 1)
          jcox= icon1[jcoxx];
        else
          jcox= icon2[jcoxx];
      
        if ( abs(jcox) != i )  {
          if ( jcox == 0 )  {
            nec_exception* nex = new nec_exception("SBF - SEGMENT CONNECTION ERROR FOR SEGMENT ");
            nex->append(i);
            throw nex;
          }
          else
            continue;
        }  
      } /* if ( jcox != i ) */
      else
      if ( jcox == is)
        *bb = - *bb;
    
      if ( iend == 1)
        break;
    } /* if ( jcox != 0 ) */
  
    pm = - pp;
    pp=0.;
    njun1= local_jsno;
  
    jcox= icon2[ix];
    if ( jcox > PCHCON)
      jcox= i;
  
    jend=1;
    iend=1;
    sig = -1.;
  
  } /* do */
  while( jcox != 0 );
  
  njun2= local_jsno- njun1;
  d= pi()* segment_length[ix];
  sdh= sin( d);
  cdh= cos( d);
  sd=2.* sdh* cdh;
  cd= cdh* cdh- sdh* sdh;
  
  if ( d <= 0.015)  {
    omc=4.* d* d;
    omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
  }
  else
    omc=1.- cd;
  
  ap=1./( log(1./( pi()* segment_radius[ix])) -.577215664);
  aj= ap;
  
  if ( njun1 == 0)  {
    if ( njun2 == 0)  {
      *aa = -1.;
      qp= pi()* segment_radius[ix];
      xxi= qp* qp;
      xxi= qp*(1.-.5* xxi)/(1.- xxi);
      *cc=1./( cdh- xxi* sdh);
        return;
    }
  
    qp= pi()* segment_radius[ix];
    xxi= qp* qp;
    xxi= qp*(1.-.5* xxi)/(1.- xxi);
    qp = -( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));
  
    if ( june == 1)  {
      *aa = - *aa* qp;
      *bb=  *bb* qp;
      *cc = - *cc* qp;
      if ( i != is)
        return;
    }
  
    *aa -= 1.;
    d = cd - xxi * sd;
    *bb += (sdh + ap * qp * (cdh - xxi * sdh)) / d;
    *cc += (cdh + ap * qp * (sdh + xxi * cdh)) / d;
    return;
  
  } /* if ( njun1 == 0) */
  
  if ( njun2 == 0)  {
    qm= pi()* segment_radius[ix];
    xxi= qm* qm;
    xxi= qm*(1.-.5* xxi)/(1.- xxi);
    qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));
  
    if ( june == -1)  {
      *aa= *aa* qm;
      *bb= *bb* qm;
      *cc= *cc* qm;
      if ( i != is)
        return;
    }
  
    *aa -= 1.;
    d= cd- xxi* sd;
    *bb += ( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
    *cc += ( cdh- aj* qm*( sdh+ xxi* cdh))/ d;
    return;  
  } /* if ( njun2 == 0) */
  
  qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
  qm=( ap* omc- pp* sd)/ qp;
  qp = -( aj* omc+ pm* sd)/ qp;
  
  if ( june != 0 )  {
    if ( june < 0 )  {
      *aa= *aa* qm;
      *bb= *bb* qm;
      *cc= *cc* qm;
    } else {
      *aa = - *aa* qp;
      *bb= *bb* qp;
      *cc = - *cc* qp;
    }
  
    if ( i != is)
      return;
  } /* if ( june != 0 ) */
  
  *aa -= 1.;
  *bb += ( aj* qm+ ap* qp)* sdh/ sd;
  *cc += ( aj* qm- ap* qp)* cdh/ sd;
}


/*
  get_current_coefficients computes coefficients of the:
    constant   [air, aii]
    sine    [bir, bii]
    cosine    [cir, cii]
  terms in the current interpolation functions for the current vector curx.
 */
void c_geometry::get_current_coefficients(nec_float wavelength, complex_array& curx,
      real_array& air, real_array& aii,
      real_array& bir, real_array& bii,
      real_array& cir, real_array& cii,
      complex_array& vqds, int nqds,
      int_array& iqds)
{
  static nec_complex s_CCJ(0.0,-0.01666666667);
  
  nec_float ar, ai, sh;
  nec_complex cs1, cs2;
  
  if ( n_segments != 0)  {
    for (int i = 0; i < n_segments; i++ )  {
      air[i] = 0.0;
      aii[i] = 0.0;
      bir[i] = 0.0;
      bii[i] = 0.0;
      cir[i] = 0.0;
      cii[i] = 0.0;
    }

    for (int i = 0; i < n_segments; i++ )  {
      ar= real( curx[i]);
      ai= imag( curx[i]);
      tbf( i+1, 1 );
    
      for (int jx = 0; jx < jsno; jx++ )  {
        int j = jco[jx]-1;
        air[j] += ax[jx]* ar;
        aii[j] += ax[jx]* ai;
        bir[j] += bx[jx]* ar;
        bii[j] += bx[jx]* ai;
        cir[j] += cx[jx]* ar;
        cii[j] += cx[jx]* ai;
      }
    } /* for( i = 0; i < n_segments; i++ ) */
  
    for (int is = 0; is < nqds; is++ )  {
      int i= iqds[is]-1;
      int jx= icon1[i];
      icon1[i]=0;
      tbf(i+1,0);
      icon1[i]= jx;
      sh = segment_length[i]*.5;
      
      
      nec_complex curd = s_CCJ * 
        vqds[is]/( (log(2.* sh/ segment_radius[i])-1.)*
        (bx[jsno-1]* cos(two_pi() * sh)+ cx[jsno-1]* sin(two_pi() * sh))* wavelength );
      ar = real( curd);
      ai = imag( curd);
    
      for ( jx = 0; jx < jsno; jx++ )  {
        int j = jco[jx]-1;
        air[j] += ax[jx]* ar;
        aii[j] += ax[jx]* ai;
        bir[j] += bx[jx]* ar;
        bii[j] += bx[jx]* ai;
        cir[j] += cx[jx]* ar;
        cii[j] += cx[jx]* ai;
      }
    
    } /* for( is = 0; is < nqds; is++ ) */
  
    for (int i = 0; i < n_segments; i++ )
      curx[i]= nec_complex( air[i]+cir[i], aii[i]+cii[i] );
  
  } /* if ( n_segments != 0) */
  
  if ( m == 0)
    return;
  
  /* convert surface currents from */
  /* t1,t2 components to x,y,z components */
  uint64_t jco1 = n_plus_2m;
  uint64_t jco2 = jco1 + m;
  
  for (int i = 1; i <= m; i++ )  {
    jco1 -= 2;
    jco2 -= 3;
    cs1= curx[jco1];
    cs2= curx[jco1+1];
    curx[jco2]  = cs1* t1x[m-i]+ cs2* t2x[m-i];
    curx[jco2+1]= cs1* t1y[m-i]+ cs2* t2y[m-i];
    curx[jco2+2]= cs1* t1z[m-i]+ cs2* t2z[m-i];
  }
}


void c_geometry::frequency_scale(nec_float freq_mhz)  {
  DEBUG_TRACE("frequency_scale(" << freq_mhz << ")");
  nec_float fr = (1.0e6 * freq_mhz) / em::speed_of_light();
  DEBUG_TRACE("       fr=(" << fr << ")");

  for (int i = 0; i < n_segments; i++ )  {
    x[i]= x_unscaled[i]* fr;
    y[i]= y_unscaled[i]* fr;
    z[i]= z_unscaled[i]* fr;
    segment_length[i]= si_unscaled[i]* fr;
    segment_radius[i]= bi_unscaled[i]* fr;
    if (segment_length[i] < 0.02)  {
      m_output->nec_printf( "WARNING- SEGMENT[%i] LENGTH TOO SMALL (%f)\n",i,segment_length[i]);
/*      nec_exception* nex = new nec_exception("SCALE - SEGMENT[");
      nex->append(i);
      nex->append("] LENGTH TOO SMALL (");
      nex->append(segment_length[i]);
      nex->append(") WAVELENGTHS ");
      throw nex;*/
    }
  }

  nec_float fr2 = fr*fr;
  for (int i = 0; i < m; i++ )  {
    px[i]= px_unscaled[i]* fr;
    py[i]= py_unscaled[i]* fr;
    pz[i]= pz_unscaled[i]* fr;
    pbi[i]= pbi_unscaled[i]* fr2;
  }
}


void c_geometry::fflds(nec_float rox, nec_float roy, nec_float roz,
      complex_array& scur, 
      nec_complex *in_ex, nec_complex *in_ey, nec_complex *in_ez ) {
  static nec_complex _const4(0.0,em::impedance() / 2.0); // +188.365
  
  // From FORTRAN common block 
  // EQUIVALENCE (XS,X), (YS,Y), (ZS,Z), (S,BI), (CONS,CONSX)
  nec_complex ex(cplx_00());
  nec_complex ey(cplx_00());
  nec_complex ez(cplx_00());
  
  for (int i = 0; i < m; i++ ) {
    nec_float arg = patch_angle(i,rox,roy,roz);
    nec_complex ct = cplx_exp(arg) * pbi[i];
    int k = 3*i;
    ex += scur[k]* ct;
    ey += scur[k+1]* ct;
    ez += scur[k+2]* ct;
  }
  
  nec_complex ct = rox*ex+ roy*ey+ roz*ez;
  
  *in_ex = _const4*(ct*rox - ex);
  *in_ey = _const4*(ct*roy - ey);
  *in_ez = _const4*(ct*roz - ez);
}


int c_geometry::test_ek_approximation(int seg1, int seg2) {
  nec_float segment_ratio = segment_radius[seg2] / segment_radius[seg1];
  
  nec_float xi = fabs(cab[seg1]*cab[seg2] + sab[seg1]*sab[seg2] + salp[seg1]*salp[seg2]);
  
  if ( (xi < 0.999999) || (fabs(segment_ratio-1.0) > 1.e-6))
    return 2;
  else
    return 0;
}

nec_float c_geometry::patch_angle(int patch_index, nec_float in_ax, nec_float in_ay, nec_float in_az) {
  return two_pi()*(in_ax*px[patch_index]+ in_ay*py[patch_index]+ in_az*pz[patch_index]);
}
  
  
