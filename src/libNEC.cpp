/***************************************************************************
 *   Copyright (C) 2004-2008,2015 by Tim Molteno                           *
 *   tim@molteno.net                                                       *
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
#include <string>

static std::string _err_message;

#define NEC_ERROR_HANDLE(__x)   { try { __x; } catch (nec_exception* _ex) {  _err_message = _ex->get_message(); return 1; } return 0;}
#define NEC_VOID_HANDLE(__x)    { try { __x; } catch (nec_exception* _ex) {  _err_message = _ex->get_message(); }}


nec_context* nec_create(void) {
  nec_context* ret = new nec_context();
  ret->initialize();
  
  return ret;
}
 
long nec_delete(nec_context* in_context) {
  NEC_ERROR_HANDLE(delete in_context);
}
 
long nec_benchmark(void) {
  return long(100.0*nec_context::benchmark());
}

long nec_wire(nec_context* in_context, int tag_id, int segment_count, double xw1, double yw1, double zw1,
  double xw2, double yw2, double zw2, double rad,
  double rdel, double rrad) {
  NEC_ERROR_HANDLE(in_context->wire(tag_id, segment_count, xw1, yw1, zw1, xw2, yw2, zw2, rad, rdel, rrad));
}
 
long nec_sp_card(nec_context* in_context, int ns,
    double x1, double y1, double z1,
    double x2, double y2, double z2) {
  NEC_ERROR_HANDLE(in_context->sp_card(ns, x1, y1, z1, x2, y2, z2));
}

long nec_sc_card(nec_context* in_context, int i2,
    double x3, double y3, double z3,
    double x4, double y4, double z4) {
  NEC_ERROR_HANDLE(in_context->sc_card(i2, x3, y3, z3, x4, y4, z4));
}

long nec_gm_card(nec_context* in_context,  int itsi, int nrpt,
                 double rox, double roy, double roz, double xs,
                 double ys, double zs, int its)
{
  double roxd = degrees_to_rad(rox);
  double royd = degrees_to_rad(roy);
  double rozd = degrees_to_rad(roz);
  NEC_ERROR_HANDLE(in_context->move(roxd, royd, rozd, xs, ys, zs, its, nrpt, itsi));
}

long nec_gx_card(nec_context* in_context, int i1, int i2) {
  NEC_ERROR_HANDLE(in_context->gx_card(i1,i2));
}


long nec_geometry_complete(nec_context* in_context, int gpflag) {
  NEC_ERROR_HANDLE(in_context->geometry_complete(gpflag));
}


const char* nec_error_message(void) {
  return _err_message.c_str();
}

double nec_gain(nec_context* in_context, int freq_index, int theta_index, int phi_index) {
  return in_context->get_gain(freq_index, theta_index, phi_index);
}

/* Statistics about the Gain distribution */
double nec_gain_max(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_max(freq_index);
}

double nec_gain_min(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_min(freq_index);
}

double nec_gain_mean(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_mean(freq_index);
}

double nec_gain_sd(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_sd(freq_index);
}

/********************** RHCP ********************************/
double nec_gain_rhcp_max(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_rhcp_max(freq_index);
}

double nec_gain_rhcp_min(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_rhcp_min(freq_index);
}

double nec_gain_rhcp_mean(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_rhcp_mean(freq_index);
}

double nec_gain_rhcp_sd(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_rhcp_sd(freq_index);
}

/********************** LHCP ********************************/
double nec_gain_lhcp_max(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_lhcp_max(freq_index);
}

double nec_gain_lhcp_min(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_lhcp_min(freq_index);
}

double nec_gain_lhcp_mean(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_lhcp_mean(freq_index);
}

double nec_gain_lhcp_sd(nec_context* in_context, int freq_index)  {
  return in_context->get_gain_lhcp_sd(freq_index);
}

/****************** IMPEDANCE CHARACTERISTICS *********************/

double nec_impedance_real(nec_context* in_context, int freq_index) {
  return in_context->get_impedance_real(freq_index);
}
double nec_impedance_imag(nec_context* in_context, int freq_index) {
  return in_context->get_impedance_imag(freq_index);
}

long nec_medium_parameters(nec_context* in_context, double permittivity, double permeability) {
  NEC_ERROR_HANDLE(in_context->medium_parameters(permittivity, permeability));
}


long nec_fr_card(nec_context* in_context, int in_ifrq, int in_nfrq, double in_freq_mhz, double in_del_freq) {
  NEC_ERROR_HANDLE(in_context->fr_card(in_ifrq, in_nfrq, in_freq_mhz, in_del_freq));
}


long nec_ld_card(nec_context* in_context, int ldtyp, int ldtag, int ldtagf, int ldtagt, double tmp1, double tmp2, double tmp3)
{
  NEC_ERROR_HANDLE(in_context->ld_card(ldtyp, ldtag, ldtagf, ldtagt, tmp1, tmp2, tmp3));
}


long nec_gn_card(nec_context* in_context, int iperf, int nradl, 
                 double epse, double sig, double tmp3, double tmp4, double tmp5, double tmp6) {
  NEC_ERROR_HANDLE(in_context->gn_card(iperf, nradl, epse, sig, tmp3, tmp4, tmp5, tmp6));
}

long nec_ex_card(nec_context* in_context, int extype, int i2, int i3, int i4, 
                 double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6) {
 NEC_ERROR_HANDLE(in_context->ex_card((enum excitation_type)extype, i2, i3, i4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6));
}

long nec_excitation_current(nec_context* in_context, double x, double y, double z, double a, double beta, double moment) {
 NEC_ERROR_HANDLE(in_context->ex_card(EXCITATION_CURRENT, 0, 0, 0, x, y, z, a, beta, moment));
}

long nec_excitation_voltage(nec_context* in_context, int tag, int segment, double v_real, double v_imag) {
 NEC_ERROR_HANDLE(in_context->ex_card(EXCITATION_VOLTAGE, tag, segment, 0, v_real, v_imag, 0.0, 0.0, 0.0, 0.0));
}

long nec_excitation_planewave(nec_context* in_context, int n_theta, int n_phi, 
                              double theta, double phi, double eta, double dtheta, double dphi, double pol_ratio) {
 NEC_ERROR_HANDLE(in_context->ex_card(EXCITATION_LINEAR, n_theta, n_phi, 0, theta, phi, eta, dtheta, dphi, pol_ratio));
}


long nec_tl_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
  NEC_ERROR_HANDLE(in_context->tl_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6));
}

long nec_nt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
  NEC_ERROR_HANDLE(in_context->nt_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6));
}

long nec_xq_card(nec_context* in_context, int itmp1)
{
  NEC_ERROR_HANDLE(in_context->xq_card(itmp1));
}

/* "gd" card, ground representation */
long nec_gd_card(nec_context* in_context, double tmp1, double tmp2, double tmp3, double tmp4) {
  NEC_ERROR_HANDLE(in_context->gd_card(tmp1, tmp2, tmp3, tmp4));
}

/* "rp" card, standard observation angle parameters */
long nec_rp_card(nec_context* in_context,
      int calc_mode, int n_theta, int n_phi,
      int output_format, int normalization, int D, int A,	
      double theta0, double phi0, double delta_theta, double delta_phi,
      double radial_distance, double gain_norm)
{
  NEC_ERROR_HANDLE(in_context->rp_card(calc_mode, n_theta, n_phi,
    output_format, normalization, D, A,
    theta0, phi0, delta_theta, delta_phi,
    radial_distance, gain_norm));
}

/* "pt" card, print control for current */
long nec_pt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4) {
  NEC_ERROR_HANDLE(in_context->pt_card(itmp1, itmp2, itmp3, itmp4));
}


/* "pq" card, print control for charge */
long nec_pq_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4) {
  NEC_ERROR_HANDLE(in_context->pq_card(itmp1, itmp2, itmp3, itmp4));
}

/* "kh" card, matrix integration limit */
long nec_kh_card(nec_context* in_context, double tmp1)
{ 
  NEC_ERROR_HANDLE(in_context->kh_card(tmp1));
}

long nec_ne_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
  NEC_ERROR_HANDLE(in_context->ne_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6));
}

long nec_nh_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6)
{
  NEC_ERROR_HANDLE(in_context->nh_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6));
}

/* "ek" card,  extended thin wire kernel option */
long nec_ek_card(nec_context* in_context, int itmp1)
{
  NEC_ERROR_HANDLE(in_context->set_extended_thin_wire_kernel(-1 != itmp1));
}


/* "cp" card, maximum coupling between antennas */
long nec_cp_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4)
{
  NEC_ERROR_HANDLE(in_context->cp_card(itmp1, itmp2, itmp3, itmp4));
}


/* "pl" card, plot flags 
	throws int on error.
*/
long nec_pl_card(nec_context* in_context, char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4)
{
  NEC_ERROR_HANDLE(in_context->pl_card(ploutput_filename, itmp1, itmp2, itmp3, itmp4));
}

