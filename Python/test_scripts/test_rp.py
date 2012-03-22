#instructions to compile :
#
#swig -c++ -python PyNEC.i 
#g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#g++ -shared -lstdc++ nec_context.o nec_output.o c_plot_card.o c_geometry.o misc.o nec_exception.o nec_ground.o c_ggrid.o matrix_algebra.o nec_radiation_pattern.o nec_structure_currents.o c_evlcom.o PyNEC_wrap.o -o _PyNEC.so

#example given for the "libnec" (modified in order to have several values for phi and theta)
#
#nec = nec_create();
#nec_wire(nec, 0, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 21, -0.042, 0.008, 0.017, -0.048, 0.021, -0.005, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 70, -0.048, 0.021, -0.005, 0.039, 0.032, -0.017, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 70, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 50, -0.042, 0.008, 0.017, 0.017, -0.015, 0.014, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 85, -0.027, 0.04, -0.031, 0.046, -0.01, 0.028, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 70, 0.017, -0.015, 0.014, -0.048, -0.038, -0.04, 0.001, 1.0, 1.0);
#nec_wire(nec, 0, 77, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0);
#nec_geometry_complete(nec, 0, 0);
#
#nec_gn_card(nec, -1,0,0.0, 0.0, 0.0,0.0, 0.0, 0.0);
#nec_ld_card(nec, 5,0,0,0,3.72e7,0.0,0.0);
#nec_pt_card(nec, -1, 0, 0, 0);
#nec_ex_card(nec, 1, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
#nec_fr_card(nec, 0, 2, 2400.0, 100.0);
#nec_rp_card(nec, 0, 3, 2, 0,5,0,0, 90.0, 90.0, 10.0, 10.0, 0.0, 0.0);
#gain = nec_get_maximum_gain(nec);

from PyNEC import *

#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add wires to the geometry
geo.wire(0, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0)
geo.wire(0, 21, -0.042, 0.008, 0.017, -0.048, 0.021, -0.005, 0.001, 1.0, 1.0)
geo.wire(0, 70, -0.048, 0.021, -0.005, 0.039, 0.032, -0.017, 0.001, 1.0, 1.0)
geo.wire(0, 70, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0)
geo.wire(0, 50, -0.042, 0.008, 0.017, 0.017, -0.015, 0.014, 0.001, 1.0, 1.0)
geo.wire(0, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0)
geo.wire(0, 85, -0.027, 0.04, -0.031, 0.046, -0.01, 0.028, 0.001, 1.0, 1.0)
geo.wire(0, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0)
geo.wire(0, 70, 0.017, -0.015, 0.014, -0.048, -0.038, -0.04, 0.001, 1.0, 1.0)
geo.wire(0, 77, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0)

#end of the geometry input
context.geometry_complete(0)

#add a "gn" card to specify the ground parameters
context.gn_card(-1, 0, 0, 0, 0, 0, 0, 0)

#add a "ld" card for "loading"
context.ld_card(5, 0, 0, 0, 3.72e7, 0.0, 0.0)

#add a "pt" card to ask for Control for Current on Wires to be printed
context.pt_card(0, 0, 1, 10)

#add a "ex" card to specify an excitation
context.ex_card(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

#add a "fr" card to specify the frequency 
context.fr_card(0, 2, 2400.0e6, 100.0e6)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(0, 3, 2, 0, 5, 0, 0, 90.0, 90.0, 10.0, 10.0, 1.0, 0.0)

#get the radiation_pattern
rp = context.get_radiation_pattern(0)
