/*
  NEC card dispatch visitor for the ANTLR 4 full NEC parser.

  Walks the parse tree and:
    - Calls c_geometry methods for geometry cards
    - Calls card handler functions (nec_card_parser.h) for program cards
    - Collects GW->GC taper state internally

  Field access:
    ctx->INT()   -> vector<TerminalNode*>  (multiple INT? fields)
    ctx->fnum() -> vector<FloatContext*>  (multiple float? fields)
    Each FloatContext wraps an INT or REAL terminal.

  Copyright (C) 2025  Timothy C.A. Molteno
  SPDX-License-Identifier: GPL-3.0-or-later
*/

#pragma once

#include <cmath>
#include <string>
#include <vector>
#include "NECFullBaseVisitor.h"
#include "c_geometry.h"
#include "nec_context.h"
#include "nec_exception.h"
#include "nec_card_parser.h"

// ---- Helpers -----------------------------------------------------------

// Value of the i-th INT terminal in a vector, or 0.0 if absent.
inline double ival(const std::vector<antlr4::tree::TerminalNode*>& v, size_t i) {
  return i < v.size() ? std::stod(v[i]->getText()) : 0.0;
}
// Value of a single optional INT/REAL terminal, or 0.0 if absent.
inline double sval(antlr4::tree::TerminalNode* t) {
  return t ? std::stod(t->getText()) : 0.0;
}

// Extract text from a float (INT | REAL) context at index i, or "0" if absent.
inline std::string ftext(const std::vector<NECFullParser::FnumContext*>& v, size_t i) {
  if (i >= v.size()) return "0";
  auto* t = v[i]->INT();
  if (!t) t = v[i]->REAL();
  return t ? t->getText() : "0";
}
inline double fval(const std::vector<NECFullParser::FnumContext*>& v, size_t i) {
  return std::stod(ftext(v, i));
}

// Build a nec_card from INT+float vectors for dispatch.
inline nec_card build_card(const std::string& mn,
                           const std::vector<antlr4::tree::TerminalNode*>& iv,
                           const std::vector<NECFullParser::FnumContext*>& fv) {
  nec_card c;
  c.mnemonic = mn;
  c.parameter_count = 0;
  for (size_t i = 0; i < 4; i++) {
    if (i < iv.size()) { c.i[i] = (int)std::stod(iv[i]->getText()); c.parameter_count++; }
    else {
      size_t fi = i - iv.size();
      if (fi < fv.size()) { c.i[i] = (int)std::stod(ftext(fv, fi)); c.parameter_count++; }
    }
  }
  for (size_t i = 0; i < 6; i++) {
    size_t fi = 4 + i;
    if (fi < iv.size()) { c.f[i] = std::stod(iv[fi]->getText()); c.parameter_count++; }
    else {
      size_t ffi = fi - iv.size();
      if (ffi < fv.size()) { c.f[i] = std::stod(ftext(fv, ffi)); c.parameter_count++; }
    }
  }
  return c;
}

// ---- Visitor -----------------------------------------------------------

class NecBuildVisitor : public NECFullBaseVisitor {
public:
  nec_context*  nec = nullptr;
  c_geometry*   geo = nullptr;
  int           nwire = 0;

  // Dispatch a program card via the handler table.
  void dispatch(const std::string& mn,
                const std::vector<antlr4::tree::TerminalNode*>& iv,
                const std::vector<NECFullParser::FnumContext*>& fv) {
    auto card = build_card(mn, iv, fv);
    auto* h = find_handler(mn);
    if (h) h->dispatch(*nec, card);
  }
  void dispatch_ints(const std::string& mn,
                     const std::vector<antlr4::tree::TerminalNode*>& iv) {
    dispatch(mn, iv, {});
  }
  void dispatch_single(const std::string& mn, antlr4::tree::TerminalNode* t) {
    std::vector<antlr4::tree::TerminalNode*> v;
    if (t) v.push_back(t);
    dispatch(mn, v, {});
  }

  // ---- Geometry cards --------------------------------------------------

  antlrcpp::Any visitGwCard(NECFullParser::GwCardContext* ctx) override {
    auto iv = ctx->INT(); auto fv = ctx->fnum();
    int   tag  = (int)ival(iv,0);  int   segs = (int)ival(iv,1);
    double x1  = fval(fv,0), y1 = fval(fv,1), z1 = fval(fv,2);
    double x2  = fval(fv,3), y2 = fval(fv,4), z2 = fval(fv,5);
    double rad = fval(fv,6);

    if (rad == 0.0) {
      _gw.active = true;  _gw.tag = tag;  _gw.segs = segs;
      _gw.x1 = x1;  _gw.y1 = y1;  _gw.z1 = z1;
      _gw.x2 = x2;  _gw.y2 = y2;  _gw.z2 = z2;
    } else {
      nwire++;
      geo->wire(tag, segs, x1, y1, z1, x2, y2, z2, rad, 1.0, 1.0);
    }
    return nullptr;
  }

  antlrcpp::Any visitGcCard(NECFullParser::GcCardContext* ctx) override {
    if (!_gw.active) throw nec_exception("GC without preceding tapered GW");
    _gw.active = false;
    auto fv = ctx->fnum();
    double rdel = fval(fv,0), rad1 = fval(fv,1), rad2 = fval(fv,2);
    if ((rad1 == 0) || (rad2 == 0))
      throw nec_exception("GEOMETRY DATA CARD ERROR");
    double rrad = pow(rad2 / rad1, 1.0 / (_gw.segs - 1));
    nwire++;
    geo->wire(_gw.tag,_gw.segs, _gw.x1,_gw.y1,_gw.z1,
              _gw.x2,_gw.y2,_gw.z2, rad1, rdel, rrad);
    return nullptr;
  }

  antlrcpp::Any visitGxCard(NECFullParser::GxCardContext* ctx) override {
    auto v = ctx->INT();
    int ix=(int)ival(v,0), iy=(int)ival(v,1);
    if (ix||iy) geo->reflect(ix, iy, 0, 0);
    return nullptr;
  }
  antlrcpp::Any visitGrCard(NECFullParser::GrCardContext* ctx) override {
    auto v = ctx->INT();
    geo->reflect(-1, 0, 0, (int)ival(v,0), (int)ival(v,1));
    return nullptr;
  }
  antlrcpp::Any visitGsCard(NECFullParser::GsCardContext* ctx) override {
    auto* f = ctx->fnum();
    geo->scale(f ? std::stod((f->INT() ? f->INT()->getText() : f->REAL()->getText())) : 0.0);
    return nullptr;
  }
  antlrcpp::Any visitGeCard(NECFullParser::GeCardContext* ctx) override {
    auto v = ctx->INT();
    geo->geometry_complete(nec, v.empty() ? 0 : (int)std::stod(v[0]->getText()));
    return nullptr;
  }
  antlrcpp::Any visitGmCard(NECFullParser::GmCardContext* ctx) override {
    auto iv = ctx->INT(); auto fv = ctx->fnum();
    geo->move(fval(fv,0), fval(fv,1), fval(fv,2),
              fval(fv,3), fval(fv,4), fval(fv,5),
              (int)ival(iv,0), (int)ival(iv,1), (int)ival(iv,0));
    return nullptr;
  }
  antlrcpp::Any visitSpCard(NECFullParser::SpCardContext* ctx) override {
    auto iv = ctx->INT(); auto fv = ctx->fnum(); int ns=(int)ival(iv,0);
    auto scv = ctx->scCard();
    if (!scv.empty()) {
      auto sf = scv[0]->fnum();
      geo->patch(ns,1, fval(fv,0),fval(fv,1),fval(fv,2),fval(fv,3),fval(fv,4),fval(fv,5),
                 fval(sf,0),fval(sf,1),fval(sf,2),fval(sf,3),fval(sf,4),fval(sf,5));
    } else {
      geo->patch(ns,1, fval(fv,0),fval(fv,1),fval(fv,2),fval(fv,3),fval(fv,4),fval(fv,5),
                 0,0,0,0,0,0);
    }
    return nullptr;
  }
  antlrcpp::Any visitSmCard(NECFullParser::SmCardContext* ctx) override {
    auto iv = ctx->INT(); auto fv = ctx->fnum(); int ns=(int)ival(iv,0);
    auto scv = ctx->scCard();
    if (!scv.empty()) {
      auto sf = scv[0]->fnum();
      geo->patch(ns,1, fval(fv,0),fval(fv,1),fval(fv,2),fval(fv,3),fval(fv,4),fval(fv,5),
                 fval(sf,0),fval(sf,1),fval(sf,2),fval(sf,3),fval(sf,4),fval(sf,5));
    } else {
      geo->patch(ns,1, fval(fv,0),fval(fv,1),fval(fv,2),fval(fv,3),fval(fv,4),fval(fv,5),
                 0,0,0,0,0,0);
    }
    return nullptr;
  }
  antlrcpp::Any visitGaCard(NECFullParser::GaCardContext* ctx) override {
    auto iv = ctx->INT(); auto fv = ctx->fnum(); nwire++;
    geo->arc((int)ival(iv,0),(int)ival(iv,1),fval(fv,0),fval(fv,1),fval(fv,2),fval(fv,3));
    return nullptr;
  }
  antlrcpp::Any visitGhCard(NECFullParser::GhCardContext* ctx) override {
    auto iv = ctx->INT(); auto fv = ctx->fnum(); nwire++;
    geo->helix((int)ival(iv,0),(int)ival(iv,1),fval(fv,0),fval(fv,1),
               fval(fv,2),fval(fv,3),fval(fv,4),fval(fv,5),fval(fv,6));
    return nullptr;
  }

  // ---- Program cards — dispatched via handler table --------------------

  antlrcpp::Any visitFrCard(NECFullParser::FrCardContext* ctx) override { dispatch("FR", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitLdCard(NECFullParser::LdCardContext* ctx) override { dispatch("LD", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitGnCard(NECFullParser::GnCardContext* ctx) override { dispatch("GN", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitExCard(NECFullParser::ExCardContext* ctx) override {
    // EX has a labeled first INT (i=INT) in the grammar — dispatch manually.
    // ctx->i is a Token*, not a TerminalNode*, so build the vector accordingly.
    std::vector<antlr4::tree::TerminalNode*> iv;
    // Synthetic TerminalNode wrapper isn't available; build card directly.
    nec_card card;
    card.mnemonic = "EX";
    card.i[0] = (int)std::stod(ctx->i->getText());
    card.parameter_count = 1;
    auto int_rest = ctx->INT();
    auto fv = ctx->fnum();
    for (size_t n = 0; n < 3 && n < int_rest.size(); n++)
      card.i[n+1] = (int)std::stod(int_rest[n]->getText()), card.parameter_count++;
    for (size_t n = 0; n < 6 && n < fv.size(); n++)
      card.f[n] = std::stod(ftext(fv, n)), card.parameter_count++;
    auto* h = find_handler("EX");
    if (h) h->dispatch(*nec, card);
    return nullptr;
  }
  antlrcpp::Any visitNtCard(NECFullParser::NtCardContext* ctx) override { dispatch("NT", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitTlCard(NECFullParser::TlCardContext* ctx) override { dispatch("TL", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitXqCard(NECFullParser::XqCardContext* ctx) override { dispatch_single("XQ", ctx->INT()); return nullptr; }
  antlrcpp::Any visitGdCard(NECFullParser::GdCardContext* ctx) override { dispatch("GD", {}, ctx->fnum()); return nullptr; }
  antlrcpp::Any visitRpCard(NECFullParser::RpCardContext* ctx) override { dispatch("RP", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitNxCard(NECFullParser::NxCardContext* ctx) override { dispatch_single("NX", ctx->INT()); return nullptr; }
  antlrcpp::Any visitPtCard(NECFullParser::PtCardContext* ctx) override { dispatch_ints("PT", ctx->INT()); return nullptr; }
  antlrcpp::Any visitKhCard(NECFullParser::KhCardContext* ctx) override { dispatch("KH", {}, ctx->fnum()); return nullptr; }
  antlrcpp::Any visitNeCard(NECFullParser::NeCardContext* ctx) override { dispatch("NE", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitNhCard(NECFullParser::NhCardContext* ctx) override { dispatch("NH", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitPqCard(NECFullParser::PqCardContext* ctx) override { dispatch("PQ", ctx->INT(), ctx->fnum()); return nullptr; }
  antlrcpp::Any visitEkCard(NECFullParser::EkCardContext* ctx) override { dispatch_single("EK", ctx->INT()); return nullptr; }
  antlrcpp::Any visitCpCard(NECFullParser::CpCardContext* ctx) override { dispatch_ints("CP", ctx->INT()); return nullptr; }
  antlrcpp::Any visitPlCard(NECFullParser::PlCardContext* ctx) override { dispatch_ints("PL", ctx->INT()); return nullptr; }
  antlrcpp::Any visitEnCard(NECFullParser::EnCardContext* ctx) override { dispatch_single("EN", ctx->INT()); return nullptr; }
  antlrcpp::Any visitWgCard(NECFullParser::WgCardContext* ctx) override { dispatch_ints("WG", ctx->INT()); return nullptr; }
  antlrcpp::Any visitMpCard(NECFullParser::MpCardContext* ctx) override { dispatch("MP", ctx->INT(), ctx->fnum()); return nullptr; }

private:
  struct GwState { bool active=false; int tag,segs; double x1,y1,z1,x2,y2,z2; } _gw;
};
