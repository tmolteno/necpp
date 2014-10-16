#include <iostream>
#include <signal.h>
#include <cstdlib>
#include <vector>
#include <string>

using namespace std;

#include "NECLexer.hpp"
#include "NECParser.hpp"
#include <antlr/ANTLRException.hpp>
#include <antlr/TokenStreamException.hpp>
#include <antlr/TokenStreamRecognitionException.hpp>
#include <antlr/CharStreamException.hpp>

#include "XGetopt.h"

int main( int argc, char **argv )
{
	nec_output_file s_output;

	nec_output_flags s_output_flags;
	string input_filename, output_filename;
	
	/* getopt() variables */
	extern char *optarg;
	
	if ( argc == 1 )
	{
		usage();
		exit(-1);
	}


	bool results_to_stdout = false;
	// allocate a new nec_context;
	nec_context s_context;
	
int option;
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
			nec_context::benchmark();
			exit(0);
		
		default: /* print usage and exit */
			usage();
			exit(-1);
		
		}	
	}

	ifstream infs(input_filename.c_str());
	FILE *output_fp=NULL;
	
	if ( !infs.good() )
	{
		string mesg = "nec2++: Could not open "  + input_filename;
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

//	secnds( &ex_timer );

	try
	{
		NECLexer lexer(infs);
		NECParser parser(lexer);
                parser.nec = &s_context;
 		parser.nec->set_output(s_output, s_output_flags);

		parser.startRule();
		return 0;
	}
//         catch(antlr::RecognitionException e)
//         {
//                 cout << "NEC++ Parse Error: Line:" << e.getLine() << endl;
//                 cout << " Error: " << e.toString() << endl;
//                 exit(1);
//         }
        catch(antlr::TokenStreamRecognitionException e)
        {
                cout << e.toString() << endl;
                exit(1);
        }
        catch(antlr::ANTLRException e)
        {
                cout << "NEC++ Parse Error: " << e.toString() << endl;
                exit(1);
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
}


