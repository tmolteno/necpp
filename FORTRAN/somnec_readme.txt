	The accompanying file SOMNEC.FOR contains the FORTRAN source code for
the SOMNEC program.  The SOMNEC program generates grids used by the NEC2
program for the Sommerfeld-Norton approximation to realistic grounds.
	I hand-coded the source code from the listing given in NOSC TD-116,
Volume I, Part II, pages 390 -- 433.  I changed the code given in the listing
to enable interactive user input of relative dielectric constant, ground
conductivity, and frequency in MHz.  The source code as given relied on the use
of a tape drive to accomplish communication of the results of the SOMNEC 
program to the NEC2 program (old-timers, takes you back, doesn't it?).
I elected to accomplish inter-program communication on the DEC VAX and Unix
workstations available to me by using a binary unformatted file.  Hence the
user is also prompted for a file name, and the grid arrays are written to the
file.  I also had to make changes in the NEC2 code, to tell NEC2 the name of
the SOMNEC file, to open the file, and to read in the file contents.  Users
wishing to employ the SOMNEC code with implementations of NEC2 or NEC81
expecting SOMNEC file input will have to determine the file type and format
expected and recode accordingly.  I eventually incorporated the SOMNEC 
program in NEC2 as a set of subroutines, and that is the version I am now
using.
	The file SOMNEC_OUTPUT.TXT contains the results of a run of SOMNEC
using the inputs in example 10 in Volume II of NOSC TD 116.  Note the
agreement with the results given on pages 154 -- 156 of Volume II.  The
time for execution of this case was 257 seconds on a DEC VAX 8600.
On a Silicon Graphics 4D/440 workstation, I obtained around 100 seconds
execution time for the same case, and I would expect tens of seconds when I 
port the code to a DEC AXP, as I will do soon.
	I hope this code is useful to you -- 73!

	George B. Christianson
	Princeton University Plasma Physics Laboratory

	gchristianson@pppl.gov

	(Amateur Radio NJ2P)

	October 5, 1994
