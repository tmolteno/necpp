#ifndef __MISC__
#define __MISC__

#include <stdio.h>
#include "math_util.h"

/* carriage return and line feed */
#define	CR	0x0d
#define	LF	0x0a

/* max length of a line read from input file */
#define	LINE_LEN	132

/*  usage()
 *
 *  prints usage information
 */
void usage(void);

/* Returns process time (user+system) BUT in _msec_ */
void secnds( nec_float *x);

/*  load_line()
 *
 *  Loads a line from a file, aborts on failure. Lines beginning
 *  with a '#' are ignored as comments. At the end of file EOF is
 *  returned.
 */
int load_line( char *buff, FILE *pfile );



#endif /* __MISC__ */
