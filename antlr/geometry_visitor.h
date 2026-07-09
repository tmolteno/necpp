/*
  Geometry construction visitor for ANTLR 4 NEC geometry parser.

  Walks the parse tree and calls c_geometry methods to build the
  antenna structure.  Handles GW→GC taper state internally.

  Copyright (C) 2025  Timothy C.A. Molteno
  SPDX-License-Identifier: GPL-3.0-or-later
*/

#pragma once

#include <cmath>
#include <string>
#include "NECGeometryBaseVisitor.h"
#include "c_geometry.h"
#include "nec_context.h"
#include "nec_exception.h"

class GeometryBuildVisitor : public NECGeometryBaseVisitor {
public:
  nec_context*  nec = nullptr;
  c_geometry*   geo = nullptr;
  int           nwire = 0;

  // Return the numeric value of the i-th NUMBER child, or 0.0 if absent.
  static double num(NECGeometryParser::GwCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GcCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GxCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GrCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GsCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GeCardContext* ctx, size_t) {
    auto* n = ctx->NUMBER();
    return n ? std::stod(n->getText()) : 0.0;
  }
  static double num(NECGeometryParser::GmCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::SpCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::SmCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::ScCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GaCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }
  static double num(NECGeometryParser::GhCardContext* ctx, size_t i) {
    return _num(ctx->NUMBER(), i);
  }

  antlrcpp::Any visitGwCard(NECGeometryParser::GwCardContext* ctx) override {
    int   tag  = (int)num(ctx, 0);
    int   segs = (int)num(ctx, 1);
    double x1  = num(ctx, 2), y1 = num(ctx, 3), z1 = num(ctx, 4);
    double x2  = num(ctx, 5), y2 = num(ctx, 6), z2 = num(ctx, 7);
    double rad = num(ctx, 8);

    if (rad == 0.0) {
      _gw.active = true;
      _gw.tag = tag;  _gw.segs = segs;
      _gw.x1 = x1;  _gw.y1 = y1;  _gw.z1 = z1;
      _gw.x2 = x2;  _gw.y2 = y2;  _gw.z2 = z2;
    } else {
      nwire++;
      geo->wire(tag, segs, x1, y1, z1, x2, y2, z2, rad, 1.0, 1.0);
    }
    return nullptr;
  }

  antlrcpp::Any visitGcCard(NECGeometryParser::GcCardContext* ctx) override {
    if (!_gw.active) throw nec_exception("GC without preceding tapered GW");
    _gw.active = false;

    double rdel = num(ctx, 2);
    double rad1 = num(ctx, 3);
    double rad2 = num(ctx, 4);

    if ((rad1 == 0) || (rad2 == 0))
      throw nec_exception("GEOMETRY DATA CARD ERROR");

    double rrad = pow(rad2 / rad1, 1.0 / (_gw.segs - 1));
    nwire++;
    geo->wire(_gw.tag, _gw.segs, _gw.x1, _gw.y1, _gw.z1,
              _gw.x2, _gw.y2, _gw.z2, rad1, rdel, rrad);
    return nullptr;
  }

  antlrcpp::Any visitGxCard(NECGeometryParser::GxCardContext* ctx) override {
    int ix = (int)num(ctx, 0), iy = (int)num(ctx, 1);
    if (ix || iy) geo->reflect(ix, iy, 0, 0);
    return nullptr;
  }

  antlrcpp::Any visitGrCard(NECGeometryParser::GrCardContext* ctx) override {
    geo->reflect(-1, 0, 0, (int)num(ctx, 0), (int)num(ctx, 1));
    return nullptr;
  }

  antlrcpp::Any visitGsCard(NECGeometryParser::GsCardContext* ctx) override {
    geo->scale(num(ctx, 0));
    return nullptr;
  }

  antlrcpp::Any visitGeCard(NECGeometryParser::GeCardContext* ctx) override {
    geo->geometry_complete(nec, (int)num(ctx, 0));
    nec->calc_prepare();
    return nullptr;
  }

  antlrcpp::Any visitGmCard(NECGeometryParser::GmCardContext* ctx) override {
    geo->move(num(ctx, 2), num(ctx, 3), num(ctx, 4),
              num(ctx, 5), num(ctx, 6), num(ctx, 7),
              (int)num(ctx, 0), (int)num(ctx, 1), (int)num(ctx, 0));
    return nullptr;
  }

  antlrcpp::Any visitSpCard(NECGeometryParser::SpCardContext* ctx) override {
    int ns = (int)num(ctx, 0);
    // Check for SC continuation
    auto* sc = ctx->scCard();
    if (sc) {
      geo->patch(ns, 1,
                 num(ctx, 1), num(ctx, 2), num(ctx, 3),
                 num(ctx, 4), num(ctx, 5), num(ctx, 6),
                 num(sc, 0), num(sc, 1), num(sc, 2),
                 num(sc, 3), num(sc, 4), num(sc, 5));
    } else {
      geo->patch(ns, 1,
                 num(ctx, 1), num(ctx, 2), num(ctx, 3),
                 num(ctx, 4), num(ctx, 5), num(ctx, 6),
                 0, 0, 0, 0, 0, 0);
    }
    return nullptr;
  }

  antlrcpp::Any visitSmCard(NECGeometryParser::SmCardContext* ctx) override {
    int ns = (int)num(ctx, 0);
    auto* sc = ctx->scCard();
    if (sc) {
      geo->patch(ns, 1,
                 num(ctx, 1), num(ctx, 2), num(ctx, 3),
                 num(ctx, 4), num(ctx, 5), num(ctx, 6),
                 num(sc, 0), num(sc, 1), num(sc, 2),
                 num(sc, 3), num(sc, 4), num(sc, 5));
    } else {
      geo->patch(ns, 1,
                 num(ctx, 1), num(ctx, 2), num(ctx, 3),
                 num(ctx, 4), num(ctx, 5), num(ctx, 6),
                 0, 0, 0, 0, 0, 0);
    }
    return nullptr;
  }

  antlrcpp::Any visitGaCard(NECGeometryParser::GaCardContext* ctx) override {
    nwire++;
    geo->arc((int)num(ctx,0), (int)num(ctx,1), num(ctx,2),
             num(ctx,3), num(ctx,4), num(ctx,5));
    return nullptr;
  }

  antlrcpp::Any visitGhCard(NECGeometryParser::GhCardContext* ctx) override {
    nwire++;
    geo->helix((int)num(ctx,0), (int)num(ctx,1),
               num(ctx,2), num(ctx,3), num(ctx,4), num(ctx,5),
               num(ctx,6), num(ctx,7), num(ctx,8));
    return nullptr;
  }

private:
  struct GwState {
    bool   active = false;
    int    tag, segs;
    double x1, y1, z1, x2, y2, z2;
  } _gw;

  static double _num(const std::vector<antlr4::tree::TerminalNode*>& nodes,
                     size_t i) {
    if (i >= nodes.size()) return 0.0;
    return std::stod(nodes[i]->getText());
  }
};
