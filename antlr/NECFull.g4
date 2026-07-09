/*
  ANTLR 4 Grammar for full NEC2 input files.

  Field types follow the NEC-2 FORTRAN I- and F/E-format:
    INT  — integer only (I-format: no decimal point, no exponent)
    REAL — floating   (F/E-format: decimal point and/or exponent)

  Each card's field sequence is typed per the NEC-2 specification.

  Copyright (C) 2025  Timothy C.A. Molteno (tim@physics.otago.ac.nz)
  SPDX-License-Identifier: GPL-3.0-or-later
*/

grammar NECFull;

options {
  language = Cpp;
}

// ---- Parser Rules ------------------------------------------------------

// Valid NEC file: one or more geometry cards, then GE, then one or more
// program cards.  CM/CE comments are skipped by the lexer.
necFile : geometrySection NEWLINE+ geCard NEWLINE+ programSection NEWLINE* EOF;

geometrySection : geometryCard (NEWLINE+ geometryCard)*;
programSection  : programCard (NEWLINE+ programCard)*;

geometryCard
  :  gwCard | gcCard | gxCard | grCard | gsCard | gmCard
  |  spCard | smCard | gaCard | ghCard;

programCard
  :  frCard | ldCard | gnCard | exCard | ntCard
  |  tlCard | xqCard | gdCard | rpCard | nxCard | ptCard
  |  khCard | neCard | nhCard | pqCard | ekCard | cpCard
  |  plCard | enCard | wgCard | mpCard;

// fnum: a numeric field — INT or REAL.  Used for FORTRAN F/E-format
// positions (coordinates, radii, etc.) where decimal notation is optional.
fnum : INT | REAL;

// ---- Geometry cards ----------------------------------------------------
//   Fields are required unless the NEC-2 manual marks them as <blank>-able.
//   Blanks default to 0 in the FORTRAN reader, modelled as ? in the grammar.

gwCard  : GW  INT INT fnum fnum fnum fnum fnum fnum fnum;
gcCard  : GC  INT? INT? fnum fnum fnum;
gxCard  : GX  INT INT;
grCard  : GR  INT INT;
gsCard  : GS  INT? INT? fnum;
geCard  : GE  INT INT?;
gmCard  : GM  INT INT fnum fnum fnum fnum fnum fnum fnum;
spCard  : SP  INT INT fnum fnum fnum fnum fnum fnum scCard*;
smCard  : SM  INT INT fnum fnum fnum fnum fnum fnum scCard*;
scCard  : SC  fnum fnum fnum fnum fnum fnum;
gaCard  : GA  INT INT fnum fnum fnum fnum;
ghCard  : GH  INT INT fnum fnum fnum fnum fnum fnum fnum;

// ---- Program cards -----------------------------------------------------
//   Field order follows parse_nec_card():  4 INT then 6 floats.
//   Cards with fewer fields omit trailing optional positions.

frCard  : FR  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;
ldCard  : LD  INT? INT? INT? INT? fnum? fnum? fnum?;
gnCard  : GN  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;
// EX — excitation.  Field layout depends on excitation type (first INT).
//   0: voltage source      EX 0  tag seg flag  v_real v_imag [norm]
//   1: linear wave         EX 1  tid t2  t3    f1 f2 f3 f4 [f5 f6]
//   2: right-hand circular EX 2  tid t2  t3    f1 f2 f3 f4 f5 f6
//   3: left-hand circular  EX 3  tid t2  t3    f1 f2 f3 f4 f5 f6
//   4: current source      EX 4  tid t2  t3    f1 f2 f3 f4 f5 f6
//   5: voltage discrete    EX 5  tag seg flag  v_real v_imag [norm]
exCard
  :  EX  i=INT
     ( {std::stoi($i.text) == 0 || std::stoi($i.text) == 5}?
       INT INT INT fnum fnum fnum?
     | {std::stoi($i.text) == 1}?
       INT INT INT fnum fnum fnum fnum fnum? fnum?
     | {std::stoi($i.text) >= 2 && std::stoi($i.text) <= 4}?
       INT INT INT fnum fnum fnum fnum fnum fnum
     )
  ;
ntCard  : NT  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;
tlCard  : TL  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;
xqCard  : XQ  INT?;
gdCard  : GD  fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum?;
rpCard  : RP  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum?;
nxCard  : NX  INT?;
ptCard  : PT  INT? INT? INT? INT?;
khCard  : KH  fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum?;
neCard  : NE  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum?;
nhCard  : NH  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum? fnum? fnum?;
pqCard  : PQ  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;
ekCard  : EK  INT?;
cpCard  : CP  INT? INT? INT? INT? INT? INT? INT? INT? INT?;
plCard  : PL  INT? INT? INT? INT? INT? INT? INT?;
enCard  : EN  INT?;
wgCard  : WG  INT? INT? INT? INT?;
mpCard  : MP  INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;

// ---- Lexer Rules ------------------------------------------------------

// Geometry mnemonics
GW : ('G'|'g')('W'|'w');  GC : ('G'|'g')('C'|'c');
GX : ('G'|'g')('X'|'x');  GR : ('G'|'g')('R'|'r');
GS : ('G'|'g')('S'|'s');  GE : ('G'|'g')('E'|'e');
GM : ('G'|'g')('M'|'m');  GA : ('G'|'g')('A'|'a');
GH : ('G'|'g')('H'|'h');
SP : ('S'|'s')('P'|'p');  SM : ('S'|'s')('M'|'m');
SC : ('S'|'s')('C'|'c');

// Program mnemonics
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

// INT: FORTRAN I-format — optional sign, digits, no decimal, no exponent.
// Must be listed before REAL so that "5" matches INT rather than REAL.
INT
  : ('+'|'-')? DIGIT+
  ;

// REAL: FORTRAN F/E-format — must contain decimal point and/or exponent.
REAL
  : ('+'|'-')?
    ( DIGIT+ '.' DIGIT* (('E'|'e') ('+'|'-')? DIGIT+)?
    | '.' DIGIT+ (('E'|'e') ('+'|'-')? DIGIT+)?
    | DIGIT+ ('E'|'e') ('+'|'-')? DIGIT+
    )
  ;

fragment DIGIT : [0-9];

NEWLINE : '\r'? '\n';
WS      : [ \t,]+ -> skip;

// CM/CE comment lines — skip entire line
CMT : ('C'|'c') ('M'|'m'|'E'|'e') ~[\r\n]* NEWLINE -> skip;
