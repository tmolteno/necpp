/* Translation to C++ by Tim Molteno

	Based on the C port by N. Kyriazis
 	Including some pieces from additional work by
		Jeroen Vreeken <pe1rxq@amsat.org>
	
	Fixed a few bugs in the process.
	
	Using the std vector library in preparation
	for moving to ATLAS for the matrix and vector
	operations.

		Debian Build Instructions
		
	apt-get install atlas3-base atlas3-headers atlas3-base-dev
	apt-get install refblas3-dev lapack3-dev lapack3-doc
	
	For more information on using LAPACK for doing 
	efficient computation, see
		http://seehuhn.de/comp/linear.html
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
#include <vector>
#include <string>

using namespace std;

#include "nec_context.h"

#ifndef _WIN32
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
	catch (nec_exception* nex)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
		s_output.line(nex->get_message().c_str());
		exit(1);
	}
	catch(...)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
 		s_output.line(" Unknown exception");
		exit(1);
	}
}

#include "c_geometry.h"
void benchmark();

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
	catch (nec_exception* nex)
	{
		cout << "NEC++ Runtime Error: " << endl;
		cout << nex->get_message().c_str() << endl;
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
	
	/* input card mnemonic list */
	/* "XT" stands for "exit", added for testing */
	#define CMD_NUM  21
	const char *atst[CMD_NUM] =
	{
		"FR", "LD", "GN", "EX", "NT", "TL", \
		"XQ", "GD", "RP", "NX", "PT", "KH", \
		"NE", "NH", "PQ", "EK", "CP", "PL", \
		"EN", "WG", "MP"
	};
	
	int itmp3, itmp2, itmp4;
	
	int ain_num;    /* ain mnemonic as a number */
	
	nec_float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	nec_float ex_timer;
	
	/* getopt() variables */
	extern char *optarg;
	int option;
	
#ifndef _WIN32
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
	
	bool results_to_stdout = false;
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
			results_to_stdout = true;
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
	if ( (input_fp = fopen(input_filename.c_str(), "r")) == NULL )
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
	if ( (output_fp = fopen(output_filename.c_str(), "w")) == NULL )
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
			throw new nec_exception("Error reading input file.");
		
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
					throw new nec_exception("Error reading input file (comments not terminated?)");
			
				/* separate card's id mnemonic */
				strncpy( ain, line_buf, 2 );
				ain[2] = '\0';
			
				/* write comment to output file */
				s_output.line(&line_buf[2]);
			}
		
			/* no "ce" card at end of comments */
			if ( strcmp(ain, "CE") != 0 )
			{
				throw new nec_exception("ERROR: INCORRECT LABEL FOR A COMMENT CARD");
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
			int itmp1;
			
			/* main input section - standard read statement - jumps */
			/* to appropriate section for specific parameter set up */
			int parameter_count = readmn(input_fp, output_fp, ain, &itmp1, &itmp2, &itmp3, &itmp4,
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
			
			/* identify card id mnemonic (except "ce" and "cm") */
			for( ain_num = 0; ain_num < CMD_NUM; ain_num++ )
				if ( strncmp( ain, atst[ain_num], 2) == 0 )
					break;
		
			/* take action according to card id mnemonic */
			switch( ain_num )
			{
			case 0: /* "fr" card, frequency parameters */
				s_context.fr_card(itmp1, itmp2, tmp1, tmp2);
				continue;
		
			case 1: /* "ld" card, loading parameters */
				s_context.ld_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3);
				continue;
		
			case 2: /* "gn" card, ground parameters under the antenna */
				s_context.gn_card(itmp1, itmp2, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; 
		
			case 3: /* "ex" card, excitation parameters */
				s_context.ex_card((enum excitation_type)itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* continue card input loop */
		
			case 4: /* "nt" card, network parameters */
				s_context.nt_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* continue card input loop */
		
			case 5: /* "tl" card, network parameters */
				if (parameter_count < 10)
				{
					nec_error_mode em(s_output);
					s_output.endl();
					s_output.line("nec2++: Missing parameters in \"TL\" card. Blank parameters should be specified as zero." );
					exit(0);
				}
				s_context.tl_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* continue card input loop */
		
			case 6: /* "xq" execute card - calc. including radiated fields */
				s_context.xq_card(itmp1);
				continue;
		
			case 7: /* "gd" card, ground representation */
				s_context.gd_card(tmp1, tmp2, tmp3, tmp4);
				continue; /* continue card input loop */
		
			case 8: /* "rp" card, standard observation angle parameters */
				{
					// pull out the XNDA parameters here...
					
					int XNDA = itmp4;
					int X = XNDA / 1000;
					int N = (XNDA / 100) % 10;
					int D = (XNDA / 10) % 10;
					int A = XNDA % 10;
					
					s_context.rp_card(itmp1, itmp2, itmp3, X, N, D, A, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				}
				continue; /* was break; followed by special code */
		
			case 9: /* "nx" card, do next job */
				next_job = true;
				continue; /* continue card input loop */
		
			case 10: /* "pt" card, print control for current */
				s_context.pt_card(itmp1, itmp2, itmp3, itmp4);
				continue; /* continue card input loop */
		
			case 11: /* "kh" card, matrix integration limit */
				s_context.kh_card(tmp1);
				continue; /* continue card input loop */
		
			case 12:  /* "ne" card, near field calculation parameters */
				s_context.ne_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue;
			
			case 13:  /* "nh" card, near field calculation parameters */
				s_context.nh_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue;
		
			case 14: /* "pq" card, print control for charge */
				s_context.pq_card(itmp1, itmp2, itmp3, itmp4);
				continue; /* continue card input loop */
		
			case 15: /* "ek" card,  extended thin wire kernel option */
				if (-1 == itmp1)
					s_context.set_extended_thin_wire_kernel(false);
				else
					s_context.set_extended_thin_wire_kernel(true);
				continue; /* continue card input loop */
		
			case 16: /* "cp" card, maximum coupling between antennas */
				s_context.cp_card(itmp1, itmp2, itmp3, itmp4);
				continue; /* continue card input loop */
		
			case 17: /* "pl" card, plot flags */
				{
					std::string ploutput_filename(input_filename);
					ploutput_filename += ".plt";
					try
					{
						s_context.pl_card(ploutput_filename.c_str(), itmp1, itmp2, itmp3, itmp4);
					}
					catch(...)
					{
						char mesg[88] = "nec2++: ";
					
						strcat( mesg, ploutput_filename.c_str() );
						perror( mesg );
						exit(-1);
					}
				}
				continue; /* continue card input loop */
		
			case 19: /* "wg" card, not supported */
				throw new nec_exception("\"WG\" card, not supported.");
		
			case 20: /* "MP" card. Material Parameters */
				s_context.medium_parameters(tmp1, tmp2);
				continue;
		
			default:
				if ( ain_num != 18 ) // EN card
				{
					throw new nec_exception("FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION.");
				}
			
				/******************************************************
				*** normal exit of nec2++ when all jobs complete ok ***
				******************************************************/
				s_context.all_jobs_completed();
				// put in here for the moment...
 				if (results_to_stdout)
					s_context.write_results(cout);
				
				/* time the process */
				secnds( &tmp1 );
				tmp1 -= ex_timer;
				fprintf( output_fp, "\n\n  TOTAL RUN TIME: %d msec", (int)tmp1 );
				
				if( input_fp != NULL )
					fclose( input_fp );
				if( output_fp != NULL )
					fclose(output_fp);
					
				return(0);
			} /* switch( ain_num ) */
		
			/*
				End of the main input section. 
			
				far_field_flag is true if last card was XQ or RP
				
				This is no longer used, but I am leaving it in here while I iron this out
				properly. simulate() should be called by the xq card and the rp card.
			*/
			ASSERT(false == ((ain_num == 6) || (ain_num == 8)));
			s_context.simulate(false);
		} /* while( ! next_job ) */
	
	} /* while(true)  */

	return(0);
} /* end of nec_main() */


/*-----------------------------------------------------------------------*/
/*!\brief Read a line and fill in the parameter values. 
\return The number of parameters read
*/
int readmn(FILE* input_fp, FILE* output_fp, 
	char *gm, int *i1, int *i2, int *i3, int *i4,
	nec_float *f1, nec_float *f2, nec_float *f3,
	nec_float *f4, nec_float *f5, nec_float *f6 )
{
	int parameter_count = 0;
	char line_buf[134];
	int line_idx;
	int n_integer_params = 4, n_float_params = 6;
	int integer_array[4] = { 0, 0, 0, 0 };
	nec_float r_array[6] = { 0., 0., 0., 0., 0., 0. };
	
	/* read a line from input file */
	int eof = load_line( line_buf, input_fp );
	
	/* get line length */
	int line_length = (int) strlen(line_buf );
	
	/* abort if card's mnemonic too short or missing */
	if ( line_length < 2 )
	{
		if (EOF == eof)
		{
			// insert an EN card if we get to an end of file
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
            strncpy( gm, "EN", 2 );
#pragma GCC diagnostic pop
            return 0;
		}
		else
		{
			fprintf( output_fp,
				"\n  COMMAND DATA CARD ERROR:"
				"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
			exit(-1);
		}
	}
	
	/* extract card's mnemonic code */
	strncpy( gm, line_buf, 2 );
	gm[2] = '\0';
	
	/* Exit if "XT" command read (for testing) */
	if ( strcmp( gm, "XT" ) == 0 )
	{
		fprintf( stderr,
			"\nnec2++: Exiting after an \"XT\" command in read_geometry_card()\n" );
		fprintf( output_fp,
			"\n\n  nec2++: Exiting after an \"XT\" command in read_geometry_card()" );
		exit(0);
	}
	
	/* Return if only mnemonic on card */
	if ( line_length == 2 )
	{
		*i1 = *i2 = *i3 = *i4 = 0;
		*f1 = *f2 = *f3 = *f4 = *f5 = *f6 = 0.0;
		return 0;
	}
	
	/* read integers from line */
	line_idx = 1;
	for (int i = 0; i < n_integer_params; i++ )
	{
		/* Find first numerical character */
		while( ((line_buf[++line_idx] <  '0')  ||
			(line_buf[  line_idx] >  '9')) &&
			(line_buf[  line_idx] != '+')  &&
			(line_buf[  line_idx] != '-') )
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= r_array[0];
			*f2= r_array[1];
			*f3= r_array[2];
			*f4= r_array[3];
			*f5= r_array[4];
			*f6= r_array[5];
			return parameter_count;
		}
		
		/* read an integer from line */
		integer_array[i] = atoi( &line_buf[line_idx] );
		parameter_count++;
		
		/* traverse numerical field to next ' ' or ',' or '\0' */
		line_idx--;
		while( (line_buf[++line_idx] != ' ') &&
			(line_buf[  line_idx] != ',') &&
			(line_buf[  line_idx] != '\0') )
		{
			/* test for non-numerical characters */
			if ( ((line_buf[line_idx] <  '0')  ||
				(line_buf[line_idx] >  '9')) &&
				(line_buf[line_idx] != '+')  &&
				(line_buf[line_idx] != '-') )
			{
				fprintf( output_fp,
				"\n  COMMAND DATA CARD \"%s\" ERROR:"
				"\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
				gm, line_buf[line_idx], (line_idx+1) );
				exit(-1);
			}		
		} /* while( (line_buff[++line_idx] ... */
		
		/* Return on end of line */
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= r_array[0];
			*f2= r_array[1];
			*f3= r_array[2];
			*f4= r_array[3];
			*f5= r_array[4];
			*f6= r_array[5];
			return parameter_count;
		}		
	} /* for( i = 0; i < n_integer_params; i++ ) */
	
	/* read nec_floats from line */
	for (int i = 0; i < n_float_params; i++ )
	{
		/* Find first numerical character */
		while( ((line_buf[++line_idx] <  '0')  ||
			(line_buf[  line_idx] >  '9')) &&
			(line_buf[  line_idx] != '+')  &&
			(line_buf[  line_idx] != '-')  &&
			(line_buf[  line_idx] != '.') )
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= r_array[0];
			*f2= r_array[1];
			*f3= r_array[2];
			*f4= r_array[3];
			*f5= r_array[4];
			*f6= r_array[5];
			return parameter_count;
		}
		
		/* read a nec_float from line */
		r_array[i] = atof( &line_buf[line_idx] );
		parameter_count++;
		
		/* traverse numerical field to next ' ' or ',' */
		line_idx--;
		while( (line_buf[++line_idx] != ' ') &&
			(line_buf[  line_idx] != ',') &&
			(line_buf[  line_idx] != '\0') )
		{
			/* test for non-numerical characters */
			if ( ((line_buf[line_idx] <  '0')  ||
				(line_buf[line_idx] >  '9')) &&
				(line_buf[line_idx] != '.')  &&
				(line_buf[line_idx] != '+')  &&
				(line_buf[line_idx] != '-')  &&
				(line_buf[line_idx] != 'E')  &&
				(line_buf[line_idx] != 'e') )
			{
				fprintf( output_fp,
				"\n  COMMAND DATA CARD \"%s\" ERROR:"
				"\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d\n",
				gm, line_buf[line_idx], (line_idx+1) );
				exit(-1);
			}		
		} /* while( (line_buff[++line_idx] ... */
		
		/* Return on end of line */
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= r_array[0];
			*f2= r_array[1];
			*f3= r_array[2];
			*f4= r_array[3];
			*f5= r_array[4];
			*f6= r_array[5];
			return parameter_count;
		}		
	} /* for( i = 0; i < n_float_params; i++ ) */
	
	*i1= integer_array[0];
	*i2= integer_array[1];
	*i3= integer_array[2];
	*i4= integer_array[3];
	*f1= r_array[0];
	*f2= r_array[1];
	*f3= r_array[2];
	*f4= r_array[3];
	*f5= r_array[4];
	*f6= r_array[5];
	
	return parameter_count;
}


/*-----------------------------------------------------------------------*/


#ifndef _WIN32
static void sig_handler(int signal )
{
	switch( signal )
	{
		case SIGINT :
			fprintf(stderr, "nec2++: exiting via user interrupt" );
			exit( signal );
	
		case SIGSEGV :
			fprintf(stderr, "nec2++: segmentation fault" );
			exit( signal );
	
		case SIGFPE :
			fprintf(stderr, "nec2++: floating point exception" );
			exit( signal );
	
		case SIGABRT :
			fprintf(stderr, "nec2++: abort signal received" );
			exit( signal );
	
		case SIGTERM :
			fprintf(stderr, "nec2++: termination request received" );
			exit( signal );
	}
}
#endif

