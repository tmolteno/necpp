PyNEC is a Python module wrapped from NEC-2. 

*** What is required to use PyNEC ? ***

You must have Python installed. The versions tested are 2.3 and 2.4 (but it may work with other versions).
You must have Numarray installed. The only version tested is 1.3.3 (like above, it may work with other versions).
You must have SWIG installed (for build only). The versions tested are 1.3.25 and 1.3.27. (like above, it may work with other versions).

If you lack any of these, here are some links that might help you :

for Python : http://www.python.org/
for Numarray : http://www.stsci.edu/resources/software_hardware/numarray
for SWIG : http://www.swig.org/ 


*** How can I build it ? ***

Linux version :
---------------

First generate the makefile.
To do so move to the PyNEC directory (let call it <PyNEC>) :

$>cd PyNEC 

Then run the setup script, which will generate the makefile (you may be asked for some pieces of informations) :

$>python setup.py

Then use this makefile (which should be called "makefile.linux2") to build PyNEC :

$>make -f makefile.linux2

To install the module (root privileges are required) :

$>make -f makefile.linux2 install


The compiler used is g++ (tested version : 3.2.3.)


*** How can I test PyNEC ? ***

PyNEC is provided with some test scripts (in the test_scripts/ directory).
You can try and run them to test PyNEC installation.
To further test it, you should install NEC-2, and use the provided Nec files (in the repository test_nec_files/) as input files.
These files are the Nec counterparts of the pythons scripts.

to get NEC-2 : http://www.physics.otago.ac.nz/research/electronics/nec/

*** Remark ***
PyNEC is provided with the sources of NEC-2. But as NEC-2 is still under work, it may evolve.
Note that as long as the signature of the different functions wrapped don't get modified, you can replace
the sources of NEC-2 (in the directory <PyNEC>/wrap/src/) ; then re-install PyNEC (you shouldn't have to re-generate the
makefile).
If the signature of functions are modified, you'll have to modify the interface file which declares the function(s) concerned.
As for the python part, it won't need to be modified as long as the names of functions remain unmodified (python knows nothing about the
number and types of the arguments).
