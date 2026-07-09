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
//   Field requirements are taken from NEC-2 Part 3 (nec2prt3.pdf).
//   A field shown without '?' is REQUIRED by the manual; optional / blank
//   fields carry '?'.  Where the first integer selects a mode that changes
//   the field layout, semantic predicates gate the alternatives (EX card).
//   Blank fields that sit *between* meaningful fields are kept as positional
//   '?' placeholders so zero-filled fixed-format files still parse; trailing
//   blank fields are omitted.

// FR — Frequency.  I1=IFRQ (0 linear, 1 multiplicative); I2=NFRQ (blank -> 1);
//   I3, I4 blank; F1=FMHZ (MHz, required); F2=DELFRQ (step / multiplier).
frCard  : FR INT INT? INT? INT? fnum fnum?;

// LD — Loading.  I1=LDTYP selects the load type and the float meaning.
//   -1  short / nullify all loads (rest of card blank)
//    0  series  RLC (ohms, H, F)            — each of F1,F2,F3 "if none, blank"
//    1  parallel RLC (ohms, H, F)           — each of F1,F2,F3 "if none, blank"
//    2  series  RLC per metre               — each of F1,F2,F3 "if none, blank"
//    3  parallel RLC per metre              — each of F1,F2,F3 "if none, blank"
//    4  resistance + reactance (ohms)       — F1=R, F2=X
//    5  wire conductivity (mhos/metre)      — F1 only
ldCard
  : LD i=INT
    ( {std::stoi($i.text) == -1}?
    | {std::stoi($i.text) >= 0 && std::stoi($i.text) <= 3}?
      INT? INT? INT? fnum? fnum? fnum?
    | {std::stoi($i.text) == 4}?
      INT? INT? INT? fnum fnum
    | {std::stoi($i.text) == 5}?
      INT? INT? INT? fnum
    )
  ;

// GN — Ground parameters.  I1=IPERF selects the ground model.
//   -1  nullify ground / free space (rest of card blank)
//    0  finite ground, reflection-coefficient approximation
//    1  perfectly conducting ground (rest of card blank)
//    2  finite ground, Sommerfeld/Norton method
//   For 0/2: I2=NRADL (radial-wire count); I3,I4 blank;
//            F1=EPSE, F2=SIG required; F3..F6 = screen / 2nd-medium params.
gnCard
  : GN i=INT
    ( {std::stoi($i.text) == -1 || std::stoi($i.text) == 1}?
    | {std::stoi($i.text) == 0 || std::stoi($i.text) == 2}?
      INT? INT? INT? fnum fnum fnum? fnum? fnum? fnum?
    )
  ;

// EX — Excitation.  Field layout depends on excitation type (I1).
//   0: voltage source (applied-E)  tag seg flag  v_r v_i [norm]
//   1: linear plane wave           nTH nPH flag  th ph eta dth [dph] [pol]
//   2: right-hand circular wave    nTH nPH flag  th ph eta dth dph pol
//   3: left-hand circular wave     nTH nPH flag  th ph eta dth dph pol
//   4: elementary current source   –  –  flag    x y z alpha beta moment
//   5: voltage source (slope disc.) tag seg flag  v_r v_i [norm]
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

// NT — Two-port network (short-circuit admittance matrix).  I1,I2 = port 1
//   (tag, seg); I3,I4 = port 2.  I2 = -1 nullifies all NT/TL connections and
//   leaves the rest of the card blank.  Otherwise F1..F6 = Re/Im of
//   Y11, Y12, Y22 (mhos), all required.
ntCard
  : NT i1=INT i2=INT
    ( {std::stoi($i2.text) == -1}?
    | INT INT fnum fnum fnum fnum fnum fnum
    )
  ;

// TL — Transmission line.  Integer fields as on NT.  F1 = characteristic
//   impedance in ohms (negative -> crossed line, required); F2 = length
//   (blank -> straight-line distance); F3..F6 = Re/Im shunt admittance at
//   ends 1 and 2.
tlCard  : TL INT INT INT INT fnum fnum? fnum? fnum? fnum? fnum?;

// XQ — Execute.  I1: 0 none, 1 XZ cut, 2 YZ cut, 3 both.  Rest blank.
xqCard  : XQ INT?;

// GD — Additional (second-medium) ground parameters.  All integers blank.
//   F1=epsr2, F2=sig2, F3=clt, F4=cht (all required); F5, F6 blank.
gdCard  : GD fnum fnum fnum fnum;

// RP — Radiation pattern.  I1=calc_mode (0 normal; 1 surface wave; 2 linear
//   cliff; 3 circular cliff; 4 radial-wire screen; 5 screen+linear cliff;
//   6 screen+circular cliff).  I2=NTH, I3=NPH (blank -> 1); I4=XNDA (blank,
//   ignored when I1=1).  F1..F4 (theta0, phi0, dtheta, dphi) drive the sweep
//   and are required.  In normal mode F5=RFLD and F6=GNOR are optional; for
//   I1=1 the point is cylindrical (z, phi): F1..F4 become (z0, phi0, dz,
//   dphi) and F5=rho is REQUIRED (>~1 wavelength).
rpCard
  : RP i=INT INT? INT? INT?
    ( {std::stoi($i.text) == 1}?
      fnum fnum fnum fnum fnum fnum?
    | {std::stoi($i.text) != 1}?
      fnum fnum fnum fnum fnum? fnum?
    )
  ;

// NX — Next structure.  Rest of card blank.
nxCard  : NX;

// PT — Print control for wire currents.  I1=IPTFLG (-2 all [default if the
//   card is omitted]; -1 suppress; 0 limited; 1 receiving pattern; 2 plus
//   normalised; 3 normalised only).  I2=IPTAG, I3=IPTAGF, I4=IPTAGT select
//   the segment range.
ptCard  : PT INT INT? INT? INT?;

// KH — Interaction approximation range.  F1=RKH (distance in wavelengths).
khCard  : KH fnum;

// NE/NH — Near electric / magnetic fields.  I1=NEAR (0 rectangular,
//   1 spherical); I2,I3,I4 = point counts in the three coordinates
//   (blank -> 1).  F1,F2,F3 = first field-point position (required);
//   F4,F5,F6 = stepping increments.
neCard  : NE INT INT? INT? INT? fnum fnum fnum fnum? fnum? fnum?;
nhCard  : NH INT INT? INT? INT? fnum fnum fnum fnum? fnum? fnum?;

// PQ — Print control for charge densities.  I1=IPTFLQ (-1 suppress [default];
//   0/blank print).  I2,I3,I4 select the segment range as on PT.
pqCard  : PQ INT? INT? INT? INT?;

// EK — Extended thin-wire kernel.  I1: blank/0 use extended kernel;
//   -1 return to standard kernel.  Rest blank.
ekCard  : EK INT?;

// CP — Maximum coupling between segment pairs.  I1=TAG1, I2=SEG1, I3=TAG2,
//   I4=SEG2 (blank tag => following field is an absolute segment number).
cpCard  : CP INT? INT? INT? INT?;

// PL — Plot/output flags (not part of NEC-2 Part 3; output driver specific).
plCard  : PL INT? INT? INT? INT? INT? INT? INT?;

// EN — End of run.  No parameters.
enCard  : EN;

// WG — Write NGF file.  No parameters.
wgCard  : WG;

// MP — Medium parameters (not part of NEC-2 Part 3; nec2++/NEC-4 extension).
mpCard  : MP INT? INT? INT? INT? fnum? fnum? fnum? fnum? fnum? fnum?;

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
