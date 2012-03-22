#example2.nec (modified in order to get several excitations) :
#
#CMEXAMPLE 2. CENTER FED LINEAR ANTENNA. 
#CM           CURRENT SLOPE DISCONTINUITY SOURCE. 
#CM           1. THIN PERFECTLY CONDUCTING WIRE 
#CE           2. THIN ALUMINUM WIRE 
#GW 0 8 0. 0. -.25 0. 0. .25 .00001 
#GE 
#FR 0 3 0 0 200. 50. 
#EX 5 0 5 1 1. 0. 50.
#EX 5 0 4 1 1. 0. 50. 
#XQ
#EN

print('beginning of the test')

print('import of the module')
from PyNEC import *

print('begininning of card input')
#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add a wire to the geometry
geo.wire(0, 8, 0, 0, -.25, 0, 0, .25, .00001, 1, 1)

#end of the geometry input
context.geometry_complete(0)

#add a "fr" card to specify the frequency
context.fr_card(0, 3, 200e6, 50e6)

#add a "ex" card to specify an excitation
context.ex_card(5, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0)

#add an other "ex" card to specify a second excitation
context.ex_card(5, 0, 4, 0, 0, 1, 0, 0, 0, 0, 0)

#add a "xq" card  to force the simulation execution
context.xq_card(0)

print('end of card input\n')
print('get the first antenna_input (there are several ones, each one corresponding to one single frequency)\n')
ai=context.get_antenna_input(0)

print('get some of the available results')
print('get the array of segments numbers')
print(ai.get_segment())

print('\nget the array of complex currents')
print(ai.get_current())

print('\nget the array of powers')
print(ai.get_power())
print('\nend of the test - you can compare the results with the ones provided by NEC-2 using test_ai.nec as the input file\n')
