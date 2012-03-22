#CMAdapted from the example of the "libnec"
#CE
#GW 0 36 -0.0001 -0.0001 -0.0001 -0.0002 -0.0002 -0.0001 0.001 
#GE 
#GN 2 0 0 0 100 50 25 10 0.7 0.6
#LD 5 0 0 0 3.72e7
#FR 0 1 0 0 2400.
#PT -1
#EX 5 0 1 0 0. 0. 0. 0. 0. 0.
#RP 2 3 2 0500 90. 90. 10. 10. 0. 0.   
#EN

print('beginning of the test')

print('import of the module')
from PyNEC import *

print('beginning of card input')
#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add wires to the geometry
geo.wire(0, 36, -0.0001, -0.0001, -0.0001, -0.0002, -0.0002, -0.0001, 0.001, 1.0, 1.0)

#end of the geometry input
context.geometry_complete(0)

#add a "gn" card to specify the ground parameters
context.gn_card(2, 0, 100., 50., 25., 10., 0.7, 0.6)

#add a "ld" card for "loading"
context.ld_card(5, 0, 0, 0, 3.72e7, 0.0, 0.0)

#add a "pt" card to ask for Control for Current on Wires to be printed
context.pt_card(-1, 0, 0, 0)

#add a "ex" card to specify an excitation
context.ex_card(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

#add a "fr" card to specify the frequency 
context.fr_card(0, 2, 2400.0e6, 0)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(2, 3, 2, 0, 5, 0, 0, 90.0, 90.0, 10.0, 10.0, 0.0, 0.0)

print('end of card input\n')
print('get the radiation pattern')
rp = context.get_radiation_pattern(0)
print('get the associated ground\n')
gr = rp.get_ground()

print('get the relative dielectric constant of medium 2')
print(gr.get_relative_dielectric_constant2())

print('\nget the cliff edge distance')
print(gr.get_cliff_edge_distance())

print('\nget the cliff type')
print(gr.get_cliff_type())

print('\nend of the test - you can compare the results with the ones provided by NEC-2 using test_rp2.nec as the input file\n')
