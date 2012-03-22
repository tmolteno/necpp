from PyNEC import *

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
context.fr_card(0, 2, 2400.0e6, 100.0e6)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(2, 3, 2, 0, 5, 0, 0, 90.0, 90.0, 10.0, 10.0, 0.0, 0.0)

#get the radiation_pattern
rp = context.get_radiation_pattern(0)

#get the associated ground
gr = rp.get_ground()
