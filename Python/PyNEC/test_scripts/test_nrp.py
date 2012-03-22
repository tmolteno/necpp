#example3.nec (modified)
#
#CMEXAMPLE 3. VERTICAL HALF WAVELENGTH ANTENNA OVER GROUND 
#CM           EXTENDED THIN WIRE KERNEL USED 
#CM           1. PERFECT GROUND 
#CM           2. IMPERFECT GROUND INCLUDING GROUND WAVE AND RECEIVING 
#CE              PATTERN CALCULATIONS 
#GW 0 9 0. 0. 2. 0. 0. 7. .3 
#GE 1 
#EK 
#FR 0 1 0 0 30. 
#EX 0 0 5 0 1. 
#GN 1 
#RP 0 10 2 1301 0. 0. 10. 90. 
#GN 0 0 0 0 6. 1.000E-03 
#RP 0 10 2 1301 0. 0. 10. 90. 
#RP 1 10 1 0 1. 0. 2. 0. 1.000E+05 
#EX 1 10 1 0 0. 0. 0. 10. 
#PT 2 0 5 5 
#XQ 
#EN

print('beginning of the test')

print('import of the module')
from PyNEC import *

print('beginning of card input')
#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add a wire to the geometry
geo.wire(0, 9, 0, 0, 2, 0, 0, 7, .3, 1, 1)

#end of the geometry input
context.geometry_complete(0)

#add a "ek" card to initiate use of the extended thin-wire kernal
context.ek_card(1)

#add a "fr" card to specify the frequency 
context.fr_card(0, 1, 30e6, 0)

#add a "gn" card to specify the ground parameters
context.gn_card(0, 0, 6., 0.001, 0, 0, 0, 0)

#add a "ex" card to specify an excitation
context.ex_card(1, 10, 3, 0, 0, 0, 0, 0, 10, 20, 0)

#add a "pt" card to control the printing of currents on wire segments
context.pt_card(2, 0, 5, 5)

#add a "xq" card  to force the simulation execution
context.xq_card(0)

print('end of card input\n')
print('get the norm_rx_pattern\n')
nrp = context.get_norm_rx_pattern(0)

print('get the array of coordinates (tuples (theta, phi))')
print(nrp.get_coordinates())

print('\nget  the array of normalized receiving gains')
print(nrp.get_gain())

print('\nend of the test - you can compare the results with the ones provided by NEC-2 using test_nrp.nec as the input file\n')
