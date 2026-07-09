/*
  ANTLR 4 Grammar for NEC2 geometry cards

  Parses the geometry section of NEC2 input files.
  Supported cards: GW, GC, GX, GR, GS, GE, GM, SP, SM, SC, GA, GH

  Clean grammar — no embedded actions.  Use GeometryBuildVisitor
  (geometry_visitor.h) for geometry construction.

  Copyright (C) 2025  Timothy C.A. Molteno (tim@physics.otago.ac.nz)
  SPDX-License-Identifier: GPL-3.0-or-later
*/

grammar NECGeometry;

options {
  language = Cpp;
}

// ---- Parser Rules ------------------------------------------------------

geometryFile : (card NEWLINE)* geCard NEWLINE* EOF;

card : gwCard | gcCard | gxCard | grCard | gsCard | gmCard
     | spCard | smCard | gaCard | ghCard;

gwCard  : GW  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER?;
gcCard  : GC  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER?;
gxCard  : GX  NUMBER? NUMBER?;
grCard  : GR  NUMBER? NUMBER?;
gsCard  : GS  NUMBER? NUMBER?;
geCard  : GE  NUMBER?;
gmCard  : GM  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER?;
spCard  : SP  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? (scCard |);
smCard  : SM  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? (scCard |);
scCard  : SC  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER?;
gaCard  : GA  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER?;
ghCard  : GH  NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER? NUMBER?;

// ---- Lexer Rules ------------------------------------------------------

GW : ('G'|'g') ('W'|'w');
GC : ('G'|'g') ('C'|'c');
GX : ('G'|'g') ('X'|'x');
GR : ('G'|'g') ('R'|'r');
GS : ('G'|'g') ('S'|'s');
GE : ('G'|'g') ('E'|'e');
GM : ('G'|'g') ('M'|'m');
SP : ('S'|'s') ('P'|'p');
SM : ('S'|'s') ('M'|'m');
SC : ('S'|'s') ('C'|'c');
GA : ('G'|'g') ('A'|'a');
GH : ('G'|'g') ('H'|'h');

NUMBER
  : ('+'|'-')?
    ( DIGIT+ ('.' DIGIT*)? (('E'|'e') ('+'|'-')? DIGIT+)?
    | '.' DIGIT+ (('E'|'e') ('+'|'-')? DIGIT+)?
    )
  ;

fragment DIGIT : [0-9];

NEWLINE      : '\r'? '\n';
WS           : [ \t]+ -> skip;

// CM/CE comment lines — skip
COMMENT_LINE : ('C'|'c') ('M'|'m'|'E'|'e') ~[\r\n]* NEWLINE -> skip;
