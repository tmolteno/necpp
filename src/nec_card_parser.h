#ifndef __nec_card_parser__
#define __nec_card_parser__

/*
  Clean NEC card parser — replaces the FORTRAN-style readmn() with
  std::string-based parsing and table-driven dispatch.

  Copyright (C) 2025  Timothy C.A. Molteno

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
*/

#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <cstdlib>

#include "common.h"

/*-----------------------------------------------------------------------*/
/* Parsed NEC card                                                       */

struct nec_card {
    std::string mnemonic;
    int    i[4] = {0,0,0,0};
    double f[6] = {0.,0.,0.,0.,0.,0.};
    int    parameter_count = 0;

    bool is(const char* m) const { return mnemonic == m; }
};

/*-----------------------------------------------------------------------*/
/* Parse a single NEC card line into a nec_card struct.                 */
/* Handles: comma or space separators, negative numbers, E-notation,    */
/* inline ' comments, leading # comment lines.                          */

inline nec_card parse_nec_card(const std::string& line) {
    nec_card card;
    if (line.size() < 2)
        return card;

    /* Check for # comment lines */
    if (line[0] == '#')
        return card;

    /* Extract mnemonic */
    card.mnemonic = line.substr(0, 2);

    /* Find the ' comment delimiter and truncate */
    std::string rest = line.substr(2);
    size_t comment_pos = rest.find('\'');
    if (comment_pos != std::string::npos)
        rest = rest.substr(0, comment_pos);

    /* Parse remaining fields: replace commas with spaces for uniform parsing */
    for (auto& ch : rest)
        if (ch == ',') ch = ' ';

    std::stringstream ss(rest);
    std::string token;

    /* First 4 fields are integers */
    for (int n = 0; n < 4; n++) {
        if (!(ss >> token)) break;
        card.i[n] = std::atoi(token.c_str());
        card.parameter_count++;
    }

    /* Next 6 fields are floats */
    for (int n = 0; n < 6; n++) {
        if (!(ss >> token)) break;
        card.f[n] = std::atof(token.c_str());
        card.parameter_count++;
    }

    return card;
}

/*-----------------------------------------------------------------------*/
/* Table-driven card dispatch (replaces atst[] + giant switch)          */

class nec_context;  // forward decl

struct card_handler {
    const char* name;
    void (*dispatch)(nec_context& ctx, const nec_card& card);
};

/* Dispatch helpers — one per card type */
void handle_fr(nec_context& ctx, const nec_card& c);
void handle_ld(nec_context& ctx, const nec_card& c);
void handle_gn(nec_context& ctx, const nec_card& c);
void handle_ex(nec_context& ctx, const nec_card& c);
void handle_nt(nec_context& ctx, const nec_card& c);
void handle_tl(nec_context& ctx, const nec_card& c);
void handle_xq(nec_context& ctx, const nec_card& c);
void handle_gd(nec_context& ctx, const nec_card& c);
void handle_rp(nec_context& ctx, const nec_card& c);
void handle_nx(nec_context& ctx, const nec_card& c);
void handle_pt(nec_context& ctx, const nec_card& c);
void handle_kh(nec_context& ctx, const nec_card& c);
void handle_ne(nec_context& ctx, const nec_card& c);
void handle_nh(nec_context& ctx, const nec_card& c);
void handle_pq(nec_context& ctx, const nec_card& c);
void handle_ek(nec_context& ctx, const nec_card& c);
void handle_cp(nec_context& ctx, const nec_card& c);
void handle_pl(nec_context& ctx, const nec_card& c);
void handle_en(nec_context& ctx, const nec_card& c);
void handle_wg(nec_context& ctx, const nec_card& c);
void handle_mp(nec_context& ctx, const nec_card& c);

/* The handler table — all 21 NEC card types in one place */
static const card_handler card_handlers[] = {
    {"FR", handle_fr}, {"LD", handle_ld}, {"GN", handle_gn},
    {"EX", handle_ex}, {"NT", handle_nt}, {"TL", handle_tl},
    {"XQ", handle_xq}, {"GD", handle_gd}, {"RP", handle_rp},
    {"NX", handle_nx}, {"PT", handle_pt}, {"KH", handle_kh},
    {"NE", handle_ne}, {"NH", handle_nh}, {"PQ", handle_pq},
    {"EK", handle_ek}, {"CP", handle_cp}, {"PL", handle_pl},
    {"EN", handle_en}, {"WG", handle_wg}, {"MP", handle_mp},
};

inline const card_handler* find_handler(const std::string& mnemonic) {
    for (const auto& h : card_handlers)
        if (h.name[0] == mnemonic[0] && h.name[1] == mnemonic[1])
            return &h;
    return nullptr;
}

#endif /* __nec_card_parser__ */
