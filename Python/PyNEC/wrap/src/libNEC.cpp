/***************************************************************************
 *   Copyright (C) 2004-2005 by Tim Molteno                                  *
 *   tim@molteno.net                                                        *
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

#include "libnecpp.h"
#include "nec_context.h"
#include "nec_exception.h"

#define NEC_ERROR_HANDLE(__x)	{ try { __x; } catch (nec_exception* ex) { return 1; }}
#define NEC_VOID_HANDLE(__x)	{ try { __x; } catch (nec_exception* ex) { }}

/*! \brief Create an nec_context and initialize it.

Note: Do NOT delete or free the nec_context yourself, rather call nec_delete()
to free memory associated with the nec simulation.
*/
nec_context* nec_create()
{
	nec_context* ret = new nec_context();
	ret->initialize();
	
	return ret;
}
 
/*! \brief Delete an nec_context.
*/
long nec_delete(nec_context* in_context)
{
	NEC_ERROR_HANDLE(delete in_context);
	return 0;
}
 
long nec_benchmark()
{
	return long(100.0*nec_context::benchmark());
}

void nec_wire(nec_context* in_context, int tag_id, int segment_count, double xw1, double yw1, double zw1,
	double xw2, double yw2, double zw2, double rad,
	double rdel, double rrad)
{
	in_context->wire(tag_id, segment_count, xw1, yw1, zw1, xw2, yw2, zw2, rad, rdel, rrad);
}
 
 
void nec_geometry_complete(nec_context* in_context, int card_int_1, int card_int_2)
{
	in_context->geometry_complete(card_int_1, card_int_2);
}

double nec_get_maximum_gain(nec_context* in_context)
{
	return in_context->get_maximum_gain();
}



/**
 * FR crd
 *	@param in_context The nec_context created with nec_create()
 *	@param in_ifrq 0 is a linear range of frequencies, 1 is a log range.
 *	@param in_nfrq The number of frequencies
 *	@param in_freq_mhz The starting frequency in MHz.
 *	@param in_del_freq The frequency step (in MHz for ifrq = 0)
 */
void nec_fr_card(nec_context* in_context, int in_ifrq, int in_nfrq, double in_freq_mhz, double in_del_freq)
{
	in_context->fr_card(in_ifrq, in_nfrq, in_freq_mhz, in_del_freq);
}


void nec_ld_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3)
{
	in_context->ld_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3);
}


/* "gn" card, ground parameters under the antenna */
void nec_gn_card(nec_context* in_context, int itmp1, int itmp2, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
	in_context->gn_card(itmp1, itmp2, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

void nec_ex_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
	in_context->ex_card((enum excitation_type)itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}



void nec_tl_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
	in_context->tl_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

void nec_nt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
	in_context->nt_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

void nec_xq_card(nec_context* in_context, int itmp1)
{
	in_context->xq_card(itmp1);
}

/* "gd" card, ground representation */
void nec_gd_card(nec_context* in_context, double tmp1, double tmp2, double tmp3, double tmp4)
{
	in_context->gd_card(tmp1, tmp2, tmp3, tmp4);
}

/* "rp" card, standard observation angle parameters */
void nec_rp_card(nec_context* in_context,
	int calc_mode, int n_theta, int n_phi,
	int output_format, int normalization, int D, int A,	
	double theta0, double phi0, double delta_theta, double delta_phi,
	double radial_distance, double gain_norm)
{
	in_context->rp_card(calc_mode, n_theta, n_phi,
		output_format, normalization, D, A,
		theta0, phi0, delta_theta, delta_phi,
		radial_distance, gain_norm);
}

	/* "pt" card, print control for current */
void nec_pt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4)
{
	in_context->pt_card(itmp1, itmp2, itmp3, itmp4);
}


	/* "pq" card, print control for charge */
void nec_pq_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4)
{
	in_context->pq_card(itmp1, itmp2, itmp3, itmp4);
}



/* "kh" card, matrix integration limit */
void nec_kh_card(nec_context* in_context, double tmp1);

void nec_ne_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
	in_context->ne_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

void nec_nh_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
	in_context->nh_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
}

/* "ek" card,  extended thin wire kernel option */
void nec_ek_card(nec_context* in_context, int itmp1)
{
	in_context->set_extended_thin_wire_kernel(-1 != itmp1);
}


/* "cp" card, maximum coupling between antennas */
void nec_cp_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4)
{
	in_context->cp_card(itmp1, itmp2, itmp3, itmp4);
}


/* "pl" card, plot flags 
	throws int on error.
*/
void nec_pl_card(nec_context* in_context, char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4)
{
	in_context->pl_card(ploutput_filename, itmp1, itmp2, itmp3, itmp4);
}

