NEC2++ a C++ port of Numerical Electromagnetic Code

This is a free (GPL) NEC-2 compatable electromagnetic code. It can both read nec2 antenna description files (like the original) and also be incorporated into other projects like GUI tools and automatic antenna optimization systems. Nec2++ is developed on debian linux, but will work on a variety of other operating systems.

==Instructions for Linux==

** NOTE ** nec2++ now requires LAPACK (unless you use the --without-lapack command line switch)

See the INSTALL file. But here is the short version
	make -f Makefile.cvs
	./configure
	make
	sudo make install

Instructions for Compiling on Windows

Versions of nec2++ since 1.2.3 now compile fine with the MinGW (http://www.mingw.org/) free compiler tools. Just download the source distribution, and follow the unix installation guide (./configure and make).

Compiling with Visual Studio 7

NEC2++ has been tested with Microsoft Visual Studio 7. A project file for VC++ 7.0 is part of the source code distribution. Visual Studio 6.0 is now supported because it does not conform to modern C++ standards.

Step-by-step instructions

1) Unzip the Windows source code distribution.
2) Build the project inside the win32 subdirectory with Visual Studio 7.
