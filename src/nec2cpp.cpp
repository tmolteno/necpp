/* Translation to C++ by Tim Molteno

	Based on the C port by N. Kyriazis
	Including some pieces from additional work by
		Jeroen Vreeken <pe1rxq@amsat.org>

	Uses Eigen 3.4.0 (bundled in src/eigen3/) for linear algebra.
	Build with: make -j4  (just g++ and make needed)
*/

/* Original disclaimer that came with the FORTRAN code */

/******* Translated to the C language by N. Kyriazis  20 Aug 2003 *******/
/*									*/
/* Program NEC(input,tape5=input,output,tape11,tape12,tape13,tape14,	*/
/* tape15,tape16,tape20,tape21)						*/
/*									*/
/* Numerical Electromagnetics Code (NEC2)  developed at Lawrence	*/
/* Livermore lab., Livermore, CA.  (contact G. Burke at 415-422-8414	*/
/* for problems with the NEC code. For problems with the vax implem- 	*/
/* entation, contact J. Breakall at 415-422-8196 or E. Domning at 415 	*/
/* 422-5936) 								*/
/* file created 4/11/80. 						*/
/*									*/
/*                ***********Notice********** 				*/
/* This computer code material was prepared as an account of work 	*/
/* sponsored by the United States government.  Neither the United 	*/
/* States nor the United States Department Of Energy, nor any of 	*/
/* their employees, nor any of their contractors, subcontractors, 	*/
/* or their employees, makes any warranty, express or implied, or	*/
/* assumes any legal liability or responsibility for the accuracy, 	*/
/* completeness or usefulness of any information, apparatus, product 	*/
/* or process disclosed, or represents that its use would not infringe 	*/
/* privately-owned rights. 						*/
/*									*/
/************************************************************************/


#include "nec2cpp.h"
#include "nec_exception.h"
#include <signal.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include "nec_context.h"

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
/* Signal handler */
static void sig_handler( int signal );
#endif


/*-------------------------------------------------------------------*/

int nec_main( int argc, char **argv, nec_output_file& s_output );

/*	New main() function

	This places an exception handler around the old main loop to
	allow errors to be nicely caught!
*/

int main( int argc, char **argv )
{
	nec_output_file s_output;


	try
	{
		nec_main(argc, argv, s_output);
	}
	catch (const char* message)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
		s_output.line(message);
		exit(1);
	}
	catch (const nec_exception& nex)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
		s_output.line(nex.get_message().c_str());
		exit(1);
	}
	catch(...)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
 		s_output.line(" Unknown exception");
		exit(1);
	}
	return 0;
}

#include "c_geometry.h"
	void benchmark();

	#include "nec_card_parser.h"

void benchmark()
{
	try {
		cout << "The nec2++ benchmark." << endl;
		cout << "nec2++ version " nec_version << endl << endl;

		nec_float bench = nec_context::benchmark();
		cout << "Your computer's score is: " << bench << " NEC's" << endl;
	}
	catch (const char* message)
	{
		cout << "NEC++ Runtime Error: " << endl;
		cout << message << endl;
		exit(1);
	}
	catch (const nec_exception& nex)
	{
		cout << "NEC++ Runtime Error: " << endl;
		cout << nex.get_message().c_str() << endl;
		exit(1);
	}
}

int readmn(FILE* input_fp, FILE* output_fp, char *gm, int *i1, int *i2, int *i3, int *i4, nec_float *f1,
	nec_float *f2, nec_float *f3, nec_float *f4, nec_float *f5, nec_float *f6);


#include "XGetopt.h"

int nec_main( int argc, char **argv, nec_output_file& s_output )
{
	nec_output_flags s_output_flags;
	FILE *input_fp=NULL;
	FILE *output_fp=NULL;

	string input_filename, output_filename;

	char ain[3], line_buf[LINE_LEN+1];

	nec_float ex_timer;
	nec_float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	int itmp1, itmp2, itmp3, itmp4;

	/* getopt() variables */
	extern char *optarg;
	int option;

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
	/*** signal handler related code ***/
	/* new and old actions for sigaction() */
	struct sigaction sa_new, sa_old;


	/* initialize new actions */
	sa_new.sa_handler = sig_handler;
	sigemptyset( &sa_new.sa_mask );
	sa_new.sa_flags = 0;

	/* register function to handle signals */
	sigaction( SIGINT,  &sa_new, &sa_old );
	sigaction( SIGSEGV, &sa_new, 0 );
	sigaction( SIGFPE,  &sa_new, 0 );
	sigaction( SIGTERM, &sa_new, 0 );
	sigaction( SIGABRT, &sa_new, 0 );
#endif

	/*** command line arguments handler ***/
	if ( argc == 1 )
	{
		usage();
		exit(-1);
	}

	bool summary_to_stdout = false;
	// allocate a new nec_context;
	nec_context s_context;

	/* process command line options */
	while( (option = XGetopt(argc, argv, "i:o:hvscxgb") ) != -1 )
	{
		switch( option )
		{
		case 'i' : /* specify input file name */
			input_filename = optarg;
			break;

		case 'o' : /* specify output file name */
			output_filename = optarg;
			break;

		case 'g': /* return only the maximum gain to stdout */
			s_output_flags.set_gain_only(true);
			break;

		case 's': /* return output to stdout */
			summary_to_stdout = true;
			break;

		case 'c': /* use CSV result data */
			s_context.set_results_format(RESULT_FORMAT_CSV);
			break;

		case 'x': /* use XML result data */
			s_context.set_results_format(RESULT_FORMAT_XML);
			break;

		case 'h' : /* print usage and exit */
			usage();
			exit(0);

		case 'v' : /* print nec2++ version */
#ifdef _MSC_VER
			cout << ( "nec2++ " nec_version ) << " compiler: " << _MSC_VER << endl;
#else
			cout << ( "nec2++ " nec_version ) << (" compiler: " __VERSION__) << endl;
#endif
			exit(0);

		case 'b' : /* Run benchmark */
			benchmark();
			exit(0);

		default: /* print usage and exit */
			usage();
			exit(-1);

		}
	}

	/*** open input file ***/
	if ( input_filename == "" )
	{
		cerr << "nec2++: -i input_filename is required. Use input_filename \"-\" for stdin.\n";
		exit(-1);
	}
	else if ( input_filename == "-" )
	{
		input_fp = stdin;
	}
	else if ( (input_fp = fopen(input_filename.c_str(), "r")) == NULL )
	{
		string mesg = "nec2++: "  + input_filename;
		perror( mesg.c_str() );
		exit(-1);
	}

	/* make an output file name if not */
	/* specified by user on invocation */
	if ( output_filename  == "" )
	{
		/* strip the input file name extension if there is one */
		output_filename = input_filename.substr(0, input_filename.find(".",0)) + ".out";
	}

	/* open output file */
	if ( output_filename == "-" )
	{
		output_fp = stdout;
	}
	else if ( (output_fp = fopen(output_filename.c_str(), "w")) == NULL )
	{
		string mesg = "nec2++: "  + output_filename;

		perror( mesg.c_str() );
		exit(-1);
	}

	s_output.set_file(output_fp);

	secnds( &ex_timer );

	s_context.set_output(s_output, s_output_flags);
	s_context.initialize();

  /* main execution loop, exits at various points */
  /* depending on error conditions or end of jobs */
	while( true )
	{

		s_output.end_section();
		s_output.set_indent(31);
		s_output.line(" __________________________________________");
		s_output.line("|                                          |");
		s_output.line("| NUMERICAL ELECTROMAGNETICS CODE (nec2++) |");
		s_output.line("| Implemented in 'C++' in Double Precision |");
		s_output.line("|        Version " nec_version "        |");
		s_output.line("|__________________________________________|");

		/* read a line from input file */
		if ( load_line(line_buf, input_fp) == EOF )
			throw nec_exception("Error reading input file.");

		/* separate card's id mnemonic */
		strncpy( ain, line_buf, 2 );
		ain[2] = '\0';

		/* If its an "XT" card, exit (used for debugging) */
		if ( strcmp(ain, "XT") == 0 )
		{
			nec_error_mode em(s_output);
			s_output.end_section();
			s_output.line("nec2++: Exiting after an \"XT\" command in main()");
			exit(0);
		}

		/* if its a "cm" or "ce" card start reading comments */
		if ( (strcmp(ain, "CM") == 0) ||
			(strcmp(ain, "CE") == 0) )
		{
			s_output.end_section();
			s_output.set_indent(31);
			s_output.line("---------------- COMMENTS ----------------");
			s_output.line(&line_buf[2]);
			while( strcmp(ain, "CM") == 0 )
			{
				/* read a line from input file */
				if ( load_line(line_buf, input_fp) == EOF )
					throw nec_exception("Error reading input file (comments not terminated?)");

				/* separate card's id mnemonic */
				strncpy( ain, line_buf, 2 );
				ain[2] = '\0';

				/* write comment to output file */
				s_output.line(&line_buf[2]);
			}

			/* no "ce" card at end of comments */
			if ( strcmp(ain, "CE") != 0 )
			{
				throw nec_exception("ERROR: INCORRECT LABEL FOR A COMMENT CARD");
			}
		}
		else
		{
			rewind( input_fp );
		}

		/* initializations etc from original fortran code */
		int data_card_count=0;

		/* set up geometry data in subroutine parse_geometry */
		c_geometry* geo = s_context.get_geometry();
		geo->parse_geometry(&s_context, input_fp);

		s_context.calc_prepare();
		s_output.end_section();

		/*
			Main input section, exits at various points
			depending on error conditions or end of job.
			This is called the card input loop.
		*/

		bool next_job = false; /* start next job (next structure) flag */
		while ( ! next_job )
		{
			/* main input section - standard read statement - jumps */
			/* to appropriate section for specific parameter set up */
			readmn(input_fp, output_fp, ain, &itmp1, &itmp2, &itmp3, &itmp4,
				&tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6 );

			/* If its an "XT" card, exit */
			if ( strcmp(ain, "XT" ) == 0 )
			{
				nec_error_mode em(s_output);
				s_output.endl();
				s_output.line("nec2++: Exiting after an \"XT\" command in main()" );
				exit(0);
			}

			data_card_count++;
			fprintf( output_fp,
				"\n*****  DATA CARD N0. %3d "
				"%s %3d %5d %5d %5d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E",
				data_card_count, ain, itmp1, itmp2, itmp3, itmp4,
				tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );

			/* NX card: start next job */
			if ( strncmp(ain, "NX", 2) == 0 ) {
				next_job = true;
				continue;
			}

			/* EN card: all jobs complete */
			if ( strncmp(ain, "EN", 2) == 0 ) {
				s_context.all_jobs_completed();
				if (summary_to_stdout)
					s_context.write_results(cout);
				secnds( &tmp1 );
				tmp1 -= ex_timer;
				fprintf( output_fp, "\n\n  TOTAL RUN TIME: %d msec", static_cast<int>(tmp1) );
				if( input_fp != NULL ) fclose( input_fp );
				if( output_fp != NULL ) fclose(output_fp);
				return(0);
			}

			/* PL card needs special handling (requires input_filename) */
			if ( strncmp(ain, "PL", 2) == 0 ) {
				std::string ploutput_filename(input_filename);
				ploutput_filename += ".plt";
				try {
					s_context.pl_card(ploutput_filename.c_str(), itmp1, itmp2, itmp3, itmp4);
				} catch(...) {
					char mesg[88] = "nec2++: ";
					strcat( mesg, ploutput_filename.c_str() );
					perror( mesg );
					exit(-1);
				}
				continue;
			}

			/* Table-driven dispatch for remaining 18 cards */
			nec_card card;
			card.mnemonic = std::string(ain, 2);
			card.i[0]=itmp1; card.i[1]=itmp2; card.i[2]=itmp3; card.i[3]=itmp4;
			card.f[0]=tmp1;  card.f[1]=tmp2;  card.f[2]=tmp3;
			card.f[3]=tmp4;  card.f[4]=tmp5;  card.f[5]=tmp6;

			auto* h = find_handler(card.mnemonic);
			if (h) {
				h->dispatch(s_context, card);
				continue;
			}

			throw nec_exception("FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION.");
		} /* while( ! next_job ) */

	} /* while(true)  */

	return(0);
} /* end of nec_main() */


/*-----------------------------------------------------------------------*/
/*!\brief Read a line and fill in the parameter values (stream version).
\return The number of parameters read
*/
int readmn(std::istream& is,
	char *gm, int *i1, int *i2, int *i3, int *i4,
	nec_float *f1, nec_float *f2, nec_float *f3,
	nec_float *f4, nec_float *f5, nec_float *f6 )
{
	char line_buf[134];

	int eof = load_line( line_buf, is );
	int line_length = static_cast<int>(strlen( line_buf ));

	if ( line_length < 2 )
	{
		if (EOF == eof) {
			memcpy( gm, "EN", 2 );
			return 0;
		}
		return -1;
	}

	strncpy( gm, line_buf, 2 );
	gm[2] = '\0';

	if ( line_length == 2 ) {
		*i1 = *i2 = *i3 = *i4 = 0;
		*f1 = *f2 = *f3 = *f4 = *f5 = *f6 = 0.0;
		return 0;
	}
	nec_card card = parse_nec_card(line_buf);
	*i1 = card.i[0]; *i2 = card.i[1]; *i3 = card.i[2]; *i4 = card.i[3];
	*f1 = card.f[0]; *f2 = card.f[1]; *f3 = card.f[2];
	*f4 = card.f[3]; *f5 = card.f[4]; *f6 = card.f[5];
	return card.parameter_count;
}

/* FILE* version: delegates to stream version after reading the line. */
int readmn(FILE* input_fp, FILE* output_fp,
	char *gm, int *i1, int *i2, int *i3, int *i4,
	nec_float *f1, nec_float *f2, nec_float *f3,
	nec_float *f4, nec_float *f5, nec_float *f6 )
{
	char line_buf[134];

	int eof = load_line( line_buf, input_fp );
	int line_length = static_cast<int>(strlen( line_buf ));

	if ( line_length < 2 )
	{
		if (EOF == eof) {
			memcpy( gm, "EN", 2 );
			return 0;
		}
		fprintf( output_fp,
			"\n  COMMAND DATA CARD ERROR:"
			"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
		exit(-1);
	}

	strncpy( gm, line_buf, 2 );
	gm[2] = '\0';

	if ( strcmp( gm, "XT" ) == 0 ) {
		cerr << "\nnec2++: Exiting after an \"XT\" command in read_geometry_card()\n";
		fprintf( output_fp, "\n\n  nec2++: Exiting after an \"XT\" command in read_geometry_card()" );
		exit(0);
	}

	/* Delegate to stream version for parsing */
	std::istringstream iss(line_buf);
	return readmn(iss, gm, i1, i2, i3, i4, f1, f2, f3, f4, f5, f6);
}


/*-----------------------------------------------------------------------*/


#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
static void sig_handler(int signal )
{
	switch( signal )
	{
		case SIGINT :
			cerr << "nec2++: exiting via user interrupt";
			exit( signal );

		case SIGSEGV :
			cerr << "nec2++: segmentation fault";
			exit( signal );

		case SIGFPE :
			cerr << "nec2++: floating point exception";
			exit( signal );

		case SIGABRT :
			cerr << "nec2++: abort signal received";
			exit( signal );

		case SIGTERM :
			cerr << "nec2++: termination request received";
			exit( signal );
	}
}
#endif

/*-----------------------------------------------------------------------*/
/* Table-driven card handlers (replaces atst[] + switch dispatch)       */
/* See nec_card_parser.h for the handler table declaration.             */
/*-----------------------------------------------------------------------*/

#include "nec_card_parser.h"

void handle_fr(nec_context& ctx, const nec_card& c) {
    ctx.fr_card(c.i[0], c.i[1], c.f[0], c.f[1]);
}
void handle_ld(nec_context& ctx, const nec_card& c) {
    ctx.ld_card(c.i[0], c.i[1], c.i[2], c.i[3], c.f[0], c.f[1], c.f[2]);
}
void handle_gn(nec_context& ctx, const nec_card& c) {
    ctx.gn_card(c.i[0], c.i[1], c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_ex(nec_context& ctx, const nec_card& c) {
    ctx.ex_card((enum excitation_type)c.i[0], c.i[1], c.i[2], c.i[3],
                c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_nt(nec_context& ctx, const nec_card& c) {
    ctx.nt_card(c.i[0], c.i[1], c.i[2], c.i[3],
                c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_tl(nec_context& ctx, const nec_card& c) {
    ctx.tl_card(c.i[0], c.i[1], c.i[2], c.i[3],
                c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_xq(nec_context& ctx, const nec_card& c) {
    ctx.xq_card(c.i[0]);
}
void handle_gd(nec_context& ctx, const nec_card& c) {
    ctx.gd_card(c.f[0], c.f[1], c.f[2], c.f[3]);
}
void handle_rp(nec_context& ctx, const nec_card& c) {
    int XNDA = c.i[3];
    ctx.rp_card(c.i[0], c.i[1], c.i[2],
                XNDA / 1000, (XNDA / 100) % 10, (XNDA / 10) % 10, XNDA % 10,
                c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_nx(nec_context&, const nec_card&) {
    /* NX sets next_job flag — handled in the dispatch loop */
}
void handle_pt(nec_context& ctx, const nec_card& c) {
    ctx.pt_card(c.i[0], c.i[1], c.i[2], c.i[3]);
}
void handle_kh(nec_context& ctx, const nec_card& c) {
    ctx.kh_card(c.f[0]);
}
void handle_ne(nec_context& ctx, const nec_card& c) {
    ctx.ne_card(c.i[0], c.i[1], c.i[2], c.i[3],
                c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_nh(nec_context& ctx, const nec_card& c) {
    ctx.nh_card(c.i[0], c.i[1], c.i[2], c.i[3],
                c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5]);
}
void handle_pq(nec_context& ctx, const nec_card& c) {
    ctx.pq_card(c.i[0], c.i[1], c.i[2], c.i[3]);
}
void handle_ek(nec_context& ctx, const nec_card& c) {
    ctx.set_extended_thin_wire_kernel(c.i[0] != -1);
}
void handle_cp(nec_context& ctx, const nec_card& c) {
    ctx.cp_card(c.i[0], c.i[1], c.i[2], c.i[3]);
}
void handle_pl(nec_context& ctx, const nec_card& c) {
    /* PL requires filename — not cleanly supported via table dispatch yet.
       Falls back to the old switch-based handler in nec_main. */
    (void)ctx; (void)c;
}
void handle_en(nec_context& ctx, const nec_card&) {
    ctx.all_jobs_completed();
}
void handle_wg(nec_context&, const nec_card&) {
    throw nec_exception("\"WG\" card, not supported.");
}
void handle_mp(nec_context& ctx, const nec_card& c) {
    ctx.medium_parameters(c.f[0], c.f[1]);
}
