#instructions to compile :
#
#swig -c++ -python PyNEC.i 
#g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#g++ -shared -lstdc++ nec_context.o nec_output.o c_plot_card.o c_geometry.o misc.o nec_exception.o nec_ground.o c_ggrid.o matrix_algebra.o nec_radiation_pattern.o nec_structure_currents.o c_evlcom.o PyNEC_wrap.o -o _PyNEC.so

#example3.nec
#
#CMEXAMPLE 3. VERTICAL HALF WAVELENGTH ANTENNA OVER GROUND 
#CM           EXTENDED THIN WIRE KERNEL USED 
#CM           1. PERFECT GROUND 
#CM           2. IMPERFECT GROUND INCLUDING GROUND WAVE AND RECEIVING 
#CE              PATTERN CALCULATIONS 
#GW 0 9 0. 0. 2. 0. 0. 7. .3 
#GE 1 
#EK 
#PT 0 0 3 4
#LD 0 0 0 0 1000 1 1 
#FR 0 1 0 0 30. 
#EX 0 0 5 0 1. 
#GN 1 
#RP 0 10 2 1301 0. 0. 10. 90. 
#GN 0 0 0 0 6. 1.000E-03 
#RP 0 10 2 1301 0. 0. 10. 90. 
#RP 1 10 1 0 1. 0. 2. 0. 1.000E+05 
#EX 2 10 1 0 0. 0. 0.1 10. 0. 0.6 
#PT 2 0 5 5 
#XQ 
#EN

print('beginning of the test')
import sys

print('import of the module')
from PyNEC import *

print('beginning of card input')
#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add wires to the geometry
geo.wire(0, 9, 0, 0, 2, 0, 0, 7, .3, 1, 1)

#end of the geometry input
context.geometry_complete(0)

#add a "ek" card to initiate use of the extended thin-wire kernal
context.ek_card(1)

#add a "pt" card to control the printing of currents on wire segments
context.pt_card(0, 0, 3, 4)

#add a "ld" card for "loading"
context.ld_card(0, 0, 0, 0, 1000, 1, 1)

#add a "fr" card to specify the frequency 
context.fr_card(0, 1, 30e6, 0)

#add a "ex" card to specify an excitation
context.ex_card(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0)

#add a "gn" card to specify the ground parameters
context.gn_card(1, 0, 0, 0, 0, 0, 0, 0)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(0, 10, 2, 1, 3, 0, 1, 0, 0, 10.0, 90.0, 0.0, 0.0)

#add a "gn" card to specify the ground parameters
context.gn_card(0, 0, 6, 1.000e-3, 0, 0, 0, 0)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(0, 10, 2, 1, 3, 0, 1, 0, 0, 10.0, 90.0, 0.0, 0.0)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(1, 10, 1, 0, 0, 0, 0, 1, 0, 2, 0, 1.000e5, 0.0)

#add a "ex" card to specify an excitation
context.ex_card(2, 10, 1, 0, 0, 0, 0, 0.1, 10, 0, 0.6)

#add a "pt" card to control the printing of currents on wire segments
context.pt_card(2, 0, 5, 5)

#add a "xq" card  to force the simulation execution
context.xq_card(0)  

print('end of card input\n')
print('get the currents from the first PT card (here there are nothing but the currents in wires printed in the standard output format)\n')
sc = context.get_structure_currents(0)

print('get the output format for the currents')
print(sc.get_current_output_format())

print('\ntry and get the array of theta angles for the printing of currents')
try:
	print(sc.get_current_theta())
	print('ERROR - the output format used is the standard one !')
	sys.exit(-1)
except Warning, msg :
	print('Warning : '+msg.__str__())
	print('SUCCESS - an exception has been raised as the output format used is the standard one')
	
print('\nget the array of z-coordinate of segment centers for the printing of currents')
print(sc.get_current_segment_center_z())

print('\nget the array of complex current in wire segments')
print(sc.get_current())

print('\ntry and get the array of wire segment numbers for the printing of charge densities')
try:
	print(sc.get_q_density_segment_number())
	print('ERROR - the printing of charge densities has not been requested !')
except Warning,msg :
	print('Warning : '+msg.__str__())
	print('SUCCESS - an exception has been raised as the printing of charge densities has not been requested')

print('\ntry and get the array of patch numbers')
try:
	print(sc.get_patch_number())
	print('ERROR - there are no patches !')
except Warning,msg :
	print('Warning : '+msg.__str__())
	print('SUCCESS - an exception has been raised as there are no patches')


print('\n\nget the currents from the second PT card - but it is the third currents result as there were 2 frequencies for the first one\n')
sc2 = context.get_structure_currents(2)

print('get the output format for the currents')
print(sc2.get_current_output_format())

print('\nget the array of theta angles for the printing of currents')
print(sc2.get_current_theta())
	
print('\ntry and get the array of x-coordinate of segment centers for the printing of currents')
try:
	print(sc2.get_current_segment_center_x())
	print('ERROR - Wrong output format !')
except Warning,msg :
	print('Warning : '+msg.__str__())
	print('SUCCESS - an exception has been raised as the output format used is the one designed for a receiving pattern')

print('\nget the array of complex current in wire segments')
print(sc2.get_current())
	
print('\nend of the test - you can compare the results with the ones provided by NEC-2 using test_structure_currents.nec as the input file\n')

