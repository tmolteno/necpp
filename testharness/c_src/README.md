# README File for nec2c

## INTRODUCTION:
nec2c is a translation of the NEC2 FORTRAN source code to the C language.
The translation was performed mostly "by hand" and a lot of modifications
to the original program were introduced in order to modernize the NEC2
and to remove as many built-in limitations as possible. The attendant
SOMNEC program was also translated to C and incorporated in nec2c as a
function so that Sommerfeld ground solutions are a part of the program.

## CHANGES:
The following is a list of the more significant changes incorporated into
nec2c during translation from FORTRAN to C:

* All GO TO constructs have been removed (all 961 of them!) and "spaghetti"
code sections untangled as far as was possible to the author. Still, a lot
of the code is not as clean and straightforward as might have been.

* Obsolete memory-saving practices (such as "equivalences" of different
variables) have been eliminated and memory-sharing variables have been
separated so that they are independent.

* All fixed-size arrays used in calculations have been replaced with
buffer pointers which are allocated memory dynamically according to the
needs of the program and the complexity of each structure's geometry.
There is a two-fold advantage in this - there is virtually no limit to
the complexity of a structure (number of segments/patches etc), and there
is no wasted memory in fixed arrays. Additionally, there is no need for
data storage/swapping between memory and files and therefore functions
relating to this activity and also the NGF form of solution have been
removed from the program.

* When a Sommerfeld finite ground solution is requested, since the 
SOMNEC program has been incorporated in nec2c there is no need to store
the ground grid data in a file and read it when running nec2c. Instead,
ground grid data are calculated as needed and for each new frequency if
frequency stepping is specified.

* The factr() and solve() functions have been modified to handle the
main matrix (cm) in untransposed form so that calculations are faster.

* The parser that reads the input file allows the two characters of the
mnemonic to be in lower case if preferred. It also allows comments to be
inserted anywhere in the input file in Unix style, e.g. all lines
beginning with a '#' are ignored.

* Operationally, nec2c differs from NEC2 in not being an interactive
application. Instead, nec2c is a non-interactive command-line application
which accepts an input file name and optionally an output file name.
If this is not specified, a name for the output file is made by stripping
any extensions from the input file name and adding a ".out" extension.
Furthermore, nec2c has the potential of being incorporated in another
application (like a GUI) after suitable modifications, allowing the
creation of a stand-alone program without the need for reading files
produced separately.

* My original motive for translating NEC2 into C was to make it easier
to modify and modernize and to change obsolete functions and usage. As
a result I have edited to some extend the format of the output file to
make it more "human readable" (e.g. provided a single space between
adjacent numbers to avoid a hard-to-read "chain" of numbers joined by
- signs) etc. In my humble opinion these changes make the output file
easier to read and possibly somewhat more presentable, although this is
likely to be a problem with applications that read the output file in a
rigid manner, based on the exact output format. I apologize for this
change if it causes such problems but my intention is to eventually
modify nec2c to be used as part of a graphical application, providing
results for graphical plots directly in its buffers.

## COMPILATION:
The nec2c package is very simple at this time and compilation basically
only requires a Linux platform with development tools installed (gcc,
make and optionally gdb and "valgrind" for debugging). To compile the
source code just type "make nec2c" in the nec2c directory and if all is 
well an executable binary (nec2c) should be produced. If gdb is not
installed, remove the -g option from the line in Makefile the reads:
CC = gcc -Wall -O3 -g

These changes can also be made if debugging is not of interest, thereby
reducing the size of the binary and speeding it as well. If desired,
nec2c can be installed (to /usr/local/bin) with "make install".

There is a double precision FORTRAN source (nec2dx.f) in this package
and this can be compiled and installed by typing "make nec2dx" in the
nec2c directory. It can be run by typing nec2dx and supplying an input
and output file name and it may be used to check nec2c's results for
bugs etc.

## USAGE
nec2c is run as a non-interactive command-line application
and is invoked in the following manner:

    nec2c -i<input-file-name> [-o<output-file-name>][-hv]
          -h: print this usage information and exit.
          -v: print nec2c version number and exit.

The -i option is always needed and it specifies the name of the input
file. The -o switch is optional and it specifies the output file name.
If not used, a name for the output file is made by stripping any
extensions from the input file name and adding a ".out" extension, e.g.
nec2c -i yagi.nec will cause nec2c to read yagi.nec as the input file
and produce yagi.out as the output file.

## BUGS!!
Translating such a complex and large program from FORTRAN to C and making
so many changes along the way is very prone to bugs in the new program.
I have fixed a lot of these by using various input files that hopefully
invoke most if not all of NEC2's functions but there must still be bugs
in nec2c that will surface with some specific combinations of "cards" in
some input file. The best way to check nec2c's results is to run nec2dx
with the same input file and compare results - there should be very close
agreement between them as nec2dx is also double-precision.

## Changes in latest (January 2005) release:
I have split nec2c.c into a number of smaller files to make it easier
to work on during bug-fixing or development.

Lately (December 2004) I used the "valgrind" (http://valgrind.kde.org)
tool to check nec2c and found two significant bugs in intrp() and subph()
which I (hopefully!) have fixed. I also fixed another bug that was
found by Tim Molteno in the netwk() routine.

If you intend to use valgrind (recommended!) to test nec2c for bugs
(mainly memory allocation/access errors) then do not use performance
enhancing C flags (e.g. do not use -Ox flags etc) otherwise you will get
false error reports.

## License:
nec2c is Public Domain, same as the original FORTRAN source.
Please keep any software you write incorporating nec2c in Public Domain
or at least use an open license like GPL or BSD.

## AUTHOR:

Neoklis Kyriazis   neoklisk@cytanet.com.cy

January 27 2004
  
