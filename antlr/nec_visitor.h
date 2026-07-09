/*
  NEC card dispatch visitor for the ANTLR 4 full NEC parser.

  Walks the parse tree and:
    - Calls c_geometry methods for geometry cards
    - Calls card handler functions (nec_card_parser.h) for program cards
    - Collects GW→GC taper state internally

  Copyright (C) 2025  Timothy C.A. Molteno
  SPDX-License-Identifier: GPL-3.0-or-later
*/

#pragma once

#include <cmath>
#include <string>
#include "NECFullBaseVisitor.h"
#include "c_geometry.h"
#include "nec_context.h"
#include "nec_exception.h"
#include "nec_card_parser.h"

class NecBuildVisitor : public NECFullBaseVisitor {
public:
  nec_context*  nec = nullptr;
  c_geometry*   geo = nullptr;
  int           nwire = 0;

  // ---- Helpers ----

  static double num(const std::vector<antlr4::tree::TerminalNode*>& nodes,
                    size_t i) {
    if (i >= nodes.size()) return 0.0;
    return std::stod(nodes[i]->getText());
  }

  // Create a nec_card from a parse-tree context's NUM tokens plus a mnemonic.
  // Handles both single NUM? (TerminalNode*) and multiple NUM? (vector).
  static nec_card make_card(const std::string& mnemonic,
                            antlr4::tree::TerminalNode* single) {
    nec_card c;
    c.mnemonic = mnemonic;
    if (single) { c.i[0] = (int)std::stod(single->getText()); c.parameter_count = 1; }
    return c;
  }
  static nec_card make_card(const std::string& mnemonic,
                            const std::vector<antlr4::tree::TerminalNode*>& nums) {
    nec_card c;
    c.mnemonic = mnemonic;
    c.parameter_count = 0;
    for (size_t n = 0; n < 4 && n < nums.size(); n++)
      c.i[n] = (int)std::stod(nums[n]->getText()), c.parameter_count++;
    for (size_t n = 4; n < 10 && n < nums.size(); n++)
      c.f[n-4] = std::stod(nums[n]->getText()), c.parameter_count++;
    return c;
  }
  static nec_card make_card(const std::string& mnemonic) {
    nec_card c; c.mnemonic = mnemonic; return c;
  }

  // Dispatch a program card via the handler table.
  void dispatch(const std::string& mn,
                const std::vector<antlr4::tree::TerminalNode*>& nums) {
    auto card = make_card(mn, nums);
    auto* h = find_handler(mn);
    if (h) h->dispatch(*nec, card);
  }

  void dispatch_single(const std::string& mn,
                       antlr4::tree::TerminalNode* single) {
    auto card = make_card(mn, single);
    auto* h = find_handler(mn);
    if (h) h->dispatch(*nec, card);
  }

  // ---- Geometry cards --------------------------------------------------

  antlrcpp::Any visitGwCard(NECFullParser::GwCardContext* ctx) override {
    auto n = ctx->NUM();
    int   tag  = (int)num(n,0);  int   segs = (int)num(n,1);
    double x1  = num(n,2), y1 = num(n,3), z1 = num(n,4);
    double x2  = num(n,5), y2 = num(n,6), z2 = num(n,7);
    double rad = num(n,8);

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
    auto n = ctx->NUM();
    double rdel = num(n,2), rad1 = num(n,3), rad2 = num(n,4);
    if ((rad1 == 0) || (rad2 == 0))
      throw nec_exception("GEOMETRY DATA CARD ERROR");
    double rrad = pow(rad2 / rad1, 1.0 / (_gw.segs - 1));
    nwire++;
    geo->wire(_gw.tag,_gw.segs, _gw.x1,_gw.y1,_gw.z1,
              _gw.x2,_gw.y2,_gw.z2, rad1, rdel, rrad);
    return nullptr;
  }

  antlrcpp::Any visitGxCard(NECFullParser::GxCardContext* ctx) override {
    auto n = ctx->NUM();
    int ix=(int)num(n,0), iy=(int)num(n,1);
    if (ix||iy) geo->reflect(ix, iy, 0, 0);
    return nullptr;
  }
  antlrcpp::Any visitGrCard(NECFullParser::GrCardContext* ctx) override {
    auto n = ctx->NUM();
    geo->reflect(-1, 0, 0, (int)num(n,0), (int)num(n,1));
    return nullptr;
  }
  antlrcpp::Any visitGsCard(NECFullParser::GsCardContext* ctx) override {
    auto n = ctx->NUM();
    geo->scale(num(n,0));
    return nullptr;
  }
  antlrcpp::Any visitGeCard(NECFullParser::GeCardContext* ctx) override {
    auto* n = ctx->NUM();
    geo->geometry_complete(nec, n ? (int)std::stod(n->getText()) : 0);
    return nullptr;
  }
  antlrcpp::Any visitGmCard(NECFullParser::GmCardContext* ctx) override {
    auto n = ctx->NUM();
    geo->move(num(n,2), num(n,3), num(n,4), num(n,5), num(n,6), num(n,7),
              (int)num(n,0), (int)num(n,1), (int)num(n,0));
    return nullptr;
  }
  antlrcpp::Any visitSpCard(NECFullParser::SpCardContext* ctx) override {
    auto n = ctx->NUM();  int ns=(int)num(n,0);
    auto* sc = ctx->scCard();
    if (sc) {
      auto sn = sc->NUM();
      geo->patch(ns,1, num(n,1),num(n,2),num(n,3),num(n,4),num(n,5),num(n,6),
                 num(sn,0),num(sn,1),num(sn,2),num(sn,3),num(sn,4),num(sn,5));
    } else {
      geo->patch(ns,1, num(n,1),num(n,2),num(n,3),num(n,4),num(n,5),num(n,6),
                 0,0,0,0,0,0);
    }
    return nullptr;
  }
  antlrcpp::Any visitSmCard(NECFullParser::SmCardContext* ctx) override {
    auto n = ctx->NUM();  int ns=(int)num(n,0);
    auto* sc = ctx->scCard();
    if (sc) {
      auto sn = sc->NUM();
      geo->patch(ns,1, num(n,1),num(n,2),num(n,3),num(n,4),num(n,5),num(n,6),
                 num(sn,0),num(sn,1),num(sn,2),num(sn,3),num(sn,4),num(sn,5));
    } else {
      geo->patch(ns,1, num(n,1),num(n,2),num(n,3),num(n,4),num(n,5),num(n,6),
                 0,0,0,0,0,0);
    }
    return nullptr;
  }
  antlrcpp::Any visitGaCard(NECFullParser::GaCardContext* ctx) override {
    auto n = ctx->NUM();  nwire++;
    geo->arc((int)num(n,0),(int)num(n,1),num(n,2),num(n,3),num(n,4),num(n,5));
    return nullptr;
  }
  antlrcpp::Any visitGhCard(NECFullParser::GhCardContext* ctx) override {
    auto n = ctx->NUM();  nwire++;
    geo->helix((int)num(n,0),(int)num(n,1),num(n,2),num(n,3),
               num(n,4),num(n,5),num(n,6),num(n,7),num(n,8));
    return nullptr;
  }

  // ---- Program cards — dispatched via handler table --------------------

  antlrcpp::Any visitFrCard(NECFullParser::FrCardContext* ctx) override { dispatch("FR", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitLdCard(NECFullParser::LdCardContext* ctx) override { dispatch("LD", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitGnCard(NECFullParser::GnCardContext* ctx) override { dispatch("GN", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitExCard(NECFullParser::ExCardContext* ctx) override { dispatch("EX", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitNtCard(NECFullParser::NtCardContext* ctx) override { dispatch("NT", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitTlCard(NECFullParser::TlCardContext* ctx) override { dispatch("TL", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitXqCard(NECFullParser::XqCardContext* ctx) override { dispatch_single("XQ", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitGdCard(NECFullParser::GdCardContext* ctx) override { dispatch("GD", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitRpCard(NECFullParser::RpCardContext* ctx) override { dispatch("RP", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitNxCard(NECFullParser::NxCardContext* ctx) override { dispatch_single("NX", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitPtCard(NECFullParser::PtCardContext* ctx) override { dispatch("PT", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitKhCard(NECFullParser::KhCardContext* ctx) override { dispatch("KH", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitNeCard(NECFullParser::NeCardContext* ctx) override { dispatch("NE", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitNhCard(NECFullParser::NhCardContext* ctx) override { dispatch("NH", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitPqCard(NECFullParser::PqCardContext* ctx) override { dispatch("PQ", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitEkCard(NECFullParser::EkCardContext* ctx) override { dispatch_single("EK", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitCpCard(NECFullParser::CpCardContext* ctx) override { dispatch("CP", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitPlCard(NECFullParser::PlCardContext* ctx) override { dispatch("PL", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitEnCard(NECFullParser::EnCardContext*) override { dispatch("EN", {}); return nullptr; }
  antlrcpp::Any visitWgCard(NECFullParser::WgCardContext* ctx) override { dispatch("WG", ctx->NUM()); return nullptr; }
  antlrcpp::Any visitMpCard(NECFullParser::MpCardContext* ctx) override { dispatch("MP", ctx->NUM()); return nullptr; }

private:
  struct GwState { bool active=false; int tag,segs; double x1,y1,z1,x2,y2,z2; } _gw;
};
