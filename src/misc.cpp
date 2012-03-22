/*
 * Miscellaneous support functions for nec2++
 */

#include "misc.h"

#include <stdio.h>

using namespace std;

/*------------------------------------------------------------------------*/

/*  usage()
 *
 *  prints usage information
 */

void usage(void)
{
  fprintf( stderr,
      "usage: nec2++ [-i<input-file-name>] [-o<output-file-name>]"
      "\n       -g: print maximum gain to stdout."
      "\n       -b: Perform NEC++ Benchmark."
      "\n       -s: print results to standard output."
      "\n       -c: print results in comma-separated-value (CSV) format,"
      "\n           this options is used in conjunction with (-s) above."
      "\n       -h: print this usage information and exit."
      "\n       -v: print nec2++ version number and exit.\n");

} /* end of usage() */


/*------------------------------------------------------------------------*/
/* Returns process time (user+system) BUT in _msec_ */

#ifndef _WIN32
#include <sys/times.h>
#include <unistd.h>
void secnds( nec_float *x)
{
	struct tms buffer;
	
	times(&buffer);
	*x = 1000.0 * ( (nec_float)(buffer.tms_utime + buffer.tms_stime) ) /
		( (nec_float) sysconf(_SC_CLK_TCK) );
}
#else
#include <time.h>
void secnds( nec_float *x)
{
	nec_float c = nec_float(clock());
	*x = (1000.0 * c) / CLOCKS_PER_SEC;
}
#endif


/*------------------------------------------------------------------*/

/*  load_line()
 *
 *  loads a line from a file, aborts on failure. lines beginning
 *  with a '#' are ignored as comments. at the end of file EOF is
 *  returned.
 */

int load_line( char *buff, FILE *pfile )
{
	int
		num_chr = 0, /* number of characters read, excluding lf/cr */
		eof = 0, /* EOF flag */
		chr;     /* character read by getc */
	
	/* clear buffer at start */
	buff[0] = '\0';
	
	/* ignore commented lines, white spaces and eol/cr */
	if( (chr = fgetc(pfile)) == EOF )
		return( EOF );
	
	while( (chr == '#') ||
			(chr == ' ') ||
			(chr == CR ) ||
			(chr == LF ) )
	{
		/* go to the end of line (look for lf or cr) */
		while( (chr != CR) && (chr != LF) )
			if( (chr = fgetc(pfile)) == EOF )
				return( EOF );
	
		/* dump any cr/lf remaining */
		while( (chr == CR) || (chr == LF) )
			if( (chr = fgetc(pfile)) == EOF )
				return( EOF );
	
	} /* end of while( (chr == '#') || ... */
	
	while( num_chr < LINE_LEN )
	{
		/* if lf/cr reached before filling buffer, return */
		if( (chr == CR) || (chr == LF) )
			break;
	
		/* enter new char to buffer */
		buff[num_chr++] = chr;
	
		/* terminate buffer as a string on EOF */
		if( (chr = fgetc(pfile)) == EOF )
		{
			buff[num_chr] = '\0';
			eof = EOF;
		}
	} /* end of while( num_chr < max_chr ) */
	
	/* Capitalize first two characters (mnemonics) */
	if( (buff[0] > 0x60) && (buff[0] < 0x79) )
		buff[0] -= 0x20;
	if( (buff[1] > 0x60) && (buff[1] < 0x79) )
		buff[1] -= 0x20;
	
	/* terminate buffer as a string */
	buff[num_chr] = '\0';
	
	return( eof );
} /* end of load_line() */

