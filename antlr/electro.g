/*
  ANTLR Grammar for nec2++
  
  This file describes a grammar for a new language for antenna modeling.
  The language is more explicit than the NEC2 card-based language.
  
  Copyright (C) 2014  Timothy C.A. Molteno (tim@physics.otago.ac.nz)
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 
 The inspiration for this language is OpenSCAD.
  
    /* This is a comment */
    w0 = wire(start=[1,2,3], end=[2,3,4], r=0.01);
    w1 = wire(start=w0.end,  end=[2,3,5], r=0.01);
    helix();
    patch();

    curve(p0=w1.end, p1=[], p2=[], p3=[], r=0.01);   /* bezier curve */
    patch(start=[2,3,5], width=1.5, height=1.5);  /* rectangle */
    scale(3.0);  /* scale all dimensions by 3.0 */

    geometry {
      segment_size = 0.1; /* should be between 1% and 5% of a wavelength */
    }

    ground {
      type=PERFECT;
    }

    excitation {
      type=VOLTAGE;
      point=w0.start;
      freq=range(start=1.575GHz, end=1.675GHz, n=5);
      extended_thin_wire_kernel = false; /* default */
    }
        
    radiation_pattern {
      mode = "normal";
      theta = range(start=0.0, end=180, n=30);
      phi = range(start=0.0, end=90, n=30);
    }

    execute();
*/

/* GMSH-like geometry
   
    lc = 0.8;
    Point(1) = {10.0,-10.0,2.0,lc};
    Point(2) = {-10.0,-10.0,2.0,lc};
    Point(3) = {-10.0,10.0,2.0,lc};
    Point(4) = {10.0,10.0,2.0,lc};
    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Line Loop(5) = {1,2,3,4};
    Plane Surface(6) = {5};

    `gmsh #{$file_geo} -2 -o #{$file_msh}`
*/
header 
{
  // The statements in this block appear in all header files
  #include <iostream>
  ANTLR_USING_NAMESPACE(std)
  ANTLR_USING_NAMESPACE(antlr)
  
  #include "from_string.h"
  #include "nec_context.h"
  #include "c_geometry.h"
  #include "nec_exception.h"
}

options
{
  language="Cpp";
}

class AntennaParser extends Parser;
options
{
  k=4;
}
{
public:
  void reportError(const std::string& s) {
      cout << "reportError" << s << endl;
  }
  void reportWarning(const std::string& s) {
      cout << "reportWarning" << s << endl;
  }

  void reportError(const RecognitionException& ex) {
      cout << "Parse Error: " << ex.toString() << endl;
  }
  
public:
  nec_context* nec;
  nec_output_file s_output;
}

startRule
  :	( electroFile
          | /* nothing */
          )
          EOF
  ;

electroFile
  :
    {
      nec->initialize();
    }
    geometrySection
    analysisSection
    {
      nec->all_jobs_completed();
      cout << "Max Gain: " << nec->get_gain_max() << endl;
    }
  ;


/////////////////////////////////////////////////////////////
//
//	Geometry Section
//

geometrySection
        {       int gpflag = 0; }
	:
		GEOMETRY(gpflag=intNum) LBRACE (geometryExpression)+ RBRACE
                {
                        c_geometry* geo = nec->get_geometry();
                        geo->geometry_complete(nec, gpflag, 0);
                        nec->calc_prepare();
                        cout << "Geometry Complete" << endl;
                }
	;
	
geometryExpression
	:	wireExpression SEMI;
	;

	
wireExpression
	{	int tag = 0, ns = 0;
		double xw1, yw1, zw1, xw2, yw2, zw2, rad;
		double rDel = 1.0, rad1 = 1.0, rad2 = 0.0;
		
	}
	:	WIRE	tag=intNum ns=intNum
			xw1=realNum yw1=realNum zw1=realNum
			xw2=realNum yw2=realNum zw2=realNum
			(	rad=realNum NEWLINE
				{ rad != 0 }?
			|	NEWLINE
				GC (INT INT)? rDel=realNum rad1=realNum rad2=realNum NEWLINE
				{
					if ( (rad1 == 0) || (rad2 == 0) )
					{
						throw new nec_exception("GEOMETRY DATA CARD ERROR" );
					}
								
					rad = rad1;
					rad1 = pow( (rad2/rad1), (1.0/(ns-1.0)) );
				}
			)
		{
			c_geometry* geo = nec->get_geometry();
                        cout << "Wire " << tag << "," << ns 
                            << ": (" << xw1 << "," << yw1 << "," << zw1 << ")"
                            << " (" << xw2 << "," << yw2 << "," << zw2 << ")" 
                            << ", r=" << rad 
                            << endl;
			geo->wire(tag, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, rDel, rad1);
		}
	;	
	
/////////////////////////////////////////////////////////////
//
//      Excitation Section
//

excitationSection
        :       EXCITATION LBRACE (excitationLine)+ RBRACE
        ;
        
excitationLine
        :       frCard
        |       ekCard
        |       exCard
        |       gnCard
        |       ldCard
        |       ntCard
        |       rpCard
        |       ptCard
        |       tlCard
        |       xqCard
        |       enCard
        ;


/////////////////////////////////////////////////////////////
//
//	Analysis Section
//

analysisSection
	:	(analysisLine)+
	;
	
analysisLine
	:	frCard
	|	ekCard
	|	exCard
	|	gnCard
	|	ldCard
	|	ntCard
	|	rpCard
	|	ptCard
	|	tlCard
	|	xqCard
	|	enCard
	;
	
	
/*
Integers
     IFRQ (I1)  - Determines the type of frequency stepping:
                  0 - linear stepping.
                  1 - multiplicative stepping.
     NFRQ (12) - Number of frequency steps. If this field is blank,
                  one is assumed.
     (I3), (I4) - Blank.
Floating Point
     FMHZ (F1) - Frequency in MegaHertz.
     DELFRQ (F2)- Frequency stepping increment. If the frequency
                  stepping is linear, this quantity is added to the
                  frequency each time. If the stepping is multiplicative,
                  this is the multiplication factor.
     (F3) ... (F6) - Blank.
*/
frCard
	:
		{
			int ifrq;
			int nfrq=1;
			double fmhz = 0.0;
			double delfrq = 0.0;
		}
		FR ifrq=intNum nfrq=intNum (intNum intNum)?
		(	fmhz=realNum  (delfrq=realNum)*
		|	fmhz=realNum
		) NEWLINE
		{	nec->fr_card(ifrq, nfrq, fmhz, delfrq);	}
	;


/*
	Should we use the extended thin wire kernel or not. EK -1 cancels.
*/
ekCard
	:	{	int ekflg=0;	}
		EK (ekflg=intNum)? NEWLINE
		{
			nec->set_extended_thin_wire_kernel(-1 != ekflg);
		}
	;
	
	
exCard
	{	int extype = 0, tag = 0, m=0, exflag = 0;
		double f1 = 0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0, f6=0.0;
		enum excitation_type t;
	}
	:	EX extype=intNum	{ std::cout << "EX: "  << extype << endl; }
	(
		{extype == 0}? tag=intNum m=intNum exflag=intNum f1=realNum (f2=realNum (f3=realNum)?)?
		{ t = EXCITATION_VOLTAGE; }
	|	{extype == 1}? INT INT INT f1=realNum f2=realNum f3=realNum f4=realNum (f5=realNum f6=realNum)?
		{ t = EXCITATION_LINEAR; }
	|	{extype == 2}? INT INT INT f1=realNum f2=realNum f3=realNum f4=realNum f5=realNum f6=realNum
		{ t = EXCITATION_CIRC_RIGHT; }
	|	{extype == 3}? INT INT INT f1=realNum f2=realNum f3=realNum f4=realNum f5=realNum f6=realNum
		{ t = EXCITATION_CIRC_LEFT; }
	|	{extype == 4}? INT INT INT f1=realNum f2=realNum f3=realNum f4=realNum f5=realNum f6=realNum
		{ t = EXCITATION_CURRENT; }
	|	{extype == 5}? tag=intNum m=intNum exflag=intNum f1=realNum (f2=realNum f3=realNum)?
		{ t = EXCITATION_VOLTAGE_DISC; }
	) NEWLINE
	{
		nec->ex_card(t, tag, m, exflag, f1, f2, f3, f4, f5, f6 );
	}
	;

gnCard
	{	int ground_type = 0, rad_wire_count=0;
		double epse=0.0, sig=0.0, f3=0.0, f4=0.0, f5=0.0, f6=0.0;
	}
	:	GN ground_type=intNum 
	(
		{ground_type == -1}? 								// nullify previous ground conditions
	|	{ground_type == 1 }? (INT)*								// perfect ground
	|	{ground_type == 0 }? rad_wire_count=intNum INT INT epse=realNum sig=realNum (f3=realNum f4=realNum f5=realNum f6=realNum)?		// finite ground
	|	{ground_type == 2 }? rad_wire_count=intNum INT INT epse=realNum sig=realNum (f3=realNum f4=realNum f5=realNum f6=realNum)?		// Sommerfeld Norton ground
	) NEWLINE
	{
		nec->gn_card(ground_type, rad_wire_count, epse, sig, f3, f4, f5, f6 );
	}
	;

/*
       Generate Cylindrical Structure (GR)
  To reproduce a structure by rotating about the Z-axis to form a
   complete cylindrical array, and to set flags so that symmetry is
   utilized in the solution.
 
    (I1) - Tag number increment.
    (I2) - Total number of times that the structure is to occur in the
           cylindrical array.
 Decimal Numbers
    The decimal number fields are not used.
*/
grCard
	{	int tag_increment = 0, n_times=0;	}
	:	GR tag_increment=intNum n_times=intNum
	{	s_output.nec_printf(
				"\n  STRUCTURE ROTATED ABOUT Z-AXIS %d TIMES"
				" - LABELS INCREMENTED BY %d\n", n_times, tag_increment );
		
		c_geometry* geo = nec->get_geometry();
		geo->reflect( -1, 0, 0, tag_increment, n_times);
	}
	;

ldCard
	{	int ldtyp = 0, tag=0, m=0, n=0;
		double zlr,zli=0.0, zlc=0.0;
	}
	:	LD ldtyp=intNum tag=intNum m=intNum n=intNum zlr=realNum (zli=realNum zlc=realNum)? NEWLINE
	{
		nec->ld_card(ldtyp, tag, m, n, zlr, zli, zlc );
	}
	;

rpCard
	{	int rpflg = 0, ntheta, nphi, xnda;
		double thets, phis,dth, dph, rfld=0.0, gnor=0.0 ;
	}
	:	RP	rpflg=intNum ntheta=intNum nphi=intNum xnda=intNum
			thets=realNum phis=realNum dth=realNum dph=realNum (rfld=realNum (gnor=realNum)?)? NEWLINE
	{
		int X = xnda / 1000;
		int N = (xnda / 100) % 10;
		int D = (xnda / 10) % 10;
		int A = xnda % 10;
		nec->rp_card(rpflg, ntheta, nphi, X,N,D,A, thets, phis, dth, dph, rfld, gnor);
	}
	;

ptCard
	{	int a,b,c,d;
	}
	:	PT a=intNum b=intNum c=intNum d=intNum NEWLINE
	;

ntCard
	{	int tag1 = 0, m1=0, tag2, m2=0;
		double y11r, y11i, y12r, y12i, y22r, y22i;
	}
	:	NT tag1=intNum m1=intNum tag2=intNum m2=intNum
			y11r=realNum y11i=realNum y12r=realNum y12i=realNum y22r=realNum y22i=realNum NEWLINE
	{
		nec->nt_card(tag1, m1, tag2, m2, y11r, y11i, y12r, y12i, y22r, y22i);
	}
	;

/* TODO: Allow y11i to be blank */
tlCard
	{	int tag1 = 0, m1=0, tag2, m2=0;
		double y11r, y11i=0.0, y12r, y12i, y22r=0.0, y22i=0.0 ;
	}
	:	TL tag1=intNum m1=intNum tag2=intNum m2=intNum
			y11r=realNum y11i=realNum y12r=realNum y12i=realNum y22r=realNum y22i=realNum NEWLINE
	{
		nec->tl_card(tag1, m1, tag2, m2, y11r, y11i, y12r, y12i, y22r, y22i);
	}
	;

xqCard
	{	int xqflag = 0;	}
	:	XQ	(xqflag=intNum)? NEWLINE
	{
		nec->xq_card(xqflag);
		nec->simulate(false);
	}
	;

enCard
	:	EN	(NEWLINE | EOF)
	{
		nec->simulate(false);
	}
	;



parameterList
        {       }
        :       LBRACKET parameterAssign (COMMA parameterAssign)* RBRACKET
        ;

parameterAssign
        :       id = Value
        {       }
        ;


protected
realNum returns [double val]
	:	r:REAL
		{	val = from_string<double>(r->getText()); }
	|	i:INT
		{	val = double(from_string<int>(i->getText())); }
	;

protected
intNum returns [int val]
	:	txt:INT
		{	val = from_string<int>(txt->getText()); }
	;





/////////////////////////////////////////////////////////////////
//
//	Lexical Analyzer, breaks the file into tokens
//	to be used by the parser.
//
//

class AntennaLexer extends Lexer;
options
{
	caseSensitive=false;
	k=2;
	charVocabulary = '\3'..'\377';
	defaultErrorHandler=false;
}
	
GEOMETRY	:	"geometry"	;
RAD_PATTERN	:	"radiation_pattern"	;
SEMI	        :	";"	;
GE	:	"ge"	;
GX	:	"gx"	;

WIRE	:	"wire"	;
HELIX	:	"helix"	;

FOR	:	"for"	;
FR	:	"fr"	;
GN	:	"gn"	;
LD	:	"ld"	;
NE	:	"ne"	;
NT	:	"nt"	;
RP	:	"rp"	;
PT	:	"pt"	;
TL	:	"tl"	;
XQ	:	"xq"	;
PQ	:	"pq"	;
EN	:	"en"	;

WS  :   (	 '\003'..'\010' 
		| '\013' 
		| '\f' 
		| '\016'.. '\037' 
		| '\177'..'\377' 
		| ' ' 
		| '\t' )+
        {  _ttype = ANTLR_USE_NAMESPACE(antlr)Token::SKIP; }
    ;
	
NEWLINE
	options { paraphrase = "end of line"; }
	:	'\n'	{ newline(); }
	|	"\r\n"	{ newline(); }
 	;

Cpp_Comment
	:	"//" (~'\n')* '\n'	//match alt andet end newline
			{ _ttype = ANTLR_USE_NAMESPACE(antlr)Token::SKIP; newline(); }
	;

// Numbers

protected
REAL
	options { paraphrase = "a real number"; }
	:;

protected
INT
	options { paraphrase = "an integer"; }
	:;

protected
DIGIT
	:       '0'..'9'
	;

	
Number
	: 
	( '+' | '-' ) ?
	(
		( DIGIT )+
		(   
			'.' ( DIGIT )*	(Exponent)?
				{ $setType(REAL);	}   
			|	{ $setType(INT);	}
		)
	|	'.' ( DIGIT )+  (Exponent)?
			{ $setType(REAL);	}
	)
	;

Exponent
	:	'e' ( '+' | '-' )? ( DIGIT )+
	;
	
