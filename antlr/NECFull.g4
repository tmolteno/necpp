/*
  ANTLR 4 Grammar for full NEC2 input files.

  Supports all 21 NEC card types:
    Geometry:  GW, GC, GX, GR, GS, GE, GM, SP, SM, SC, GA, GH
    Program:   FR, LD, GN, EX, NT, TL, XQ, GD, RP, NX, PT, KH,
               NE, NH, PQ, EK, CP, PL, EN, WG, MP
    Comments:  CM, CE

  Use NecBuildVisitor (nec_visitor.h) for card dispatch.

  Copyright (C) 2025  Timothy C.A. Molteno (tim@physics.otago.ac.nz)
  SPDX-License-Identifier: GPL-3.0-or-later
*/

grammar NECFull;

options {
  language = Cpp;
}

// ---- Parser Rules ------------------------------------------------------

necFile : (card NEWLINE)* EOF;

card : geometryCard | programCard;

geometryCard
  :  gwCard | gcCard | gxCard | grCard | gsCard | gmCard
  |  spCard | smCard | gaCard | ghCard;

programCard
  :  geCard | frCard | ldCard | gnCard | exCard | ntCard
  |  tlCard | xqCard | gdCard | rpCard | nxCard | ptCard
  |  khCard | neCard | nhCard | pqCard | ekCard | cpCard
  |  plCard | enCard | wgCard | mpCard;

// CM/CE comments — skipped by lexer

gwCard  : GW  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
gcCard  : GC  NUM? NUM? NUM? NUM? NUM?;
gxCard  : GX  NUM? NUM?;
grCard  : GR  NUM? NUM?;
gsCard  : GS  NUM? NUM?;
geCard  : GE  NUM?;
gmCard  : GM  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
spCard  : SP  NUM? NUM? NUM? NUM? NUM? NUM? NUM? (scCard |);
smCard  : SM  NUM? NUM? NUM? NUM? NUM? NUM? NUM? (scCard |);
scCard  : SC  NUM? NUM? NUM? NUM? NUM? NUM?;
gaCard  : GA  NUM? NUM? NUM? NUM? NUM? NUM?;
ghCard  : GH  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;

// ---- Program cards -----------------------------------------------------

frCard  : FR  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
ldCard  : LD  NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
gnCard  : GN  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
exCard  : EX  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
ntCard  : NT  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
tlCard  : TL  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
xqCard  : XQ  NUM?;
gdCard  : GD  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
rpCard  : RP  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
nxCard  : NX  NUM?;
ptCard  : PT  NUM? NUM? NUM? NUM?;
khCard  : KH  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
neCard  : NE  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
nhCard  : NH  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
pqCard  : PQ  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
ekCard  : EK  NUM?;
cpCard  : CP  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
plCard  : PL  NUM? NUM? NUM? NUM? NUM? NUM? NUM?;
enCard  : EN;
wgCard  : WG  NUM? NUM? NUM? NUM?;
mpCard  : MP  NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM? NUM?;

// ---- Lexer Rules ------------------------------------------------------

// Geometry
GW : ('G'|'g')('W'|'w');  GC : ('G'|'g')('C'|'c');
GX : ('G'|'g')('X'|'x');  GR : ('G'|'g')('R'|'r');
GS : ('G'|'g')('S'|'s');  GE : ('G'|'g')('E'|'e');
GM : ('G'|'g')('M'|'m');  GA : ('G'|'g')('A'|'a');
GH : ('G'|'g')('H'|'h');
SP : ('S'|'s')('P'|'p');  SM : ('S'|'s')('M'|'m');
SC : ('S'|'s')('C'|'c');

// Program
FR : ('F'|'f')('R'|'r');  LD : ('L'|'l')('D'|'d');
GN : ('G'|'g')('N'|'n');  EX : ('E'|'e')('X'|'x');
NT : ('N'|'n')('T'|'t');  TL : ('T'|'t')('L'|'l');
XQ : ('X'|'x')('Q'|'q');  GD : ('G'|'g')('D'|'d');
RP : ('R'|'r')('P'|'p');  NX : ('N'|'n')('X'|'x');
PT : ('P'|'p')('T'|'t');  KH : ('K'|'k')('H'|'h');
NE : ('N'|'n')('E'|'e');  NH : ('N'|'n')('H'|'h');
PQ : ('P'|'p')('Q'|'q');  EK : ('E'|'e')('K'|'k');
CP : ('C'|'c')('P'|'p');  PL : ('P'|'p')('L'|'l');
EN : ('E'|'e')('N'|'n');  WG : ('W'|'w')('G'|'g');
MP : ('M'|'m')('P'|'p');

NUM
  : ('+'|'-')?
    ( DIGIT+ ('.' DIGIT*)? (('E'|'e') ('+'|'-')? DIGIT+)?
    | '.' DIGIT+ (('E'|'e') ('+'|'-')? DIGIT+)?
    )
  ;

fragment DIGIT : [0-9];

NEWLINE : '\r'? '\n';
WS      : [ \t]+ -> skip;

// CM/CE comment lines — skip entire line
CMT : ('C'|'c') ('M'|'m'|'E'|'e') ~[\r\n]* NEWLINE -> skip;
