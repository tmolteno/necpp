#instructions to compile :
#
#swig -c++ -python PyNEC.i 
#g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#g++ -shared -lstdc++ nec_context.o nec_output.o c_plot_card.o c_geometry.o misc.o nec_exception.o nec_ground.o c_ggrid.o matrix_algebra.o nec_radiation_pattern.o nec_structure_currents.o c_evlcom.o PyNEC_wrap.o -o _PyNEC.so

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

from PyNEC import *

#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add a wire to the geometry
geo.wire(0, 8, 0, 0, -.25, 0, 0, .25, .00001, 1, 1)

#end of the geometry input
context.geometry_complete(0)

#add a "fr" card to specify the frequency 
context.fr_card(0, 3, 200e6, 50)

#add a "ex" card to specify an excitation
context.ex_card(5, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0)

#add an other "ex" card to specify a second excitation
context.ex_card(5, 0, 4, 0, 0, 1, 0, 0, 0, 0, 0)

#add a "xq" card  to force the simulation execution
context.xq_card(0)

#get the first antenna_input (there are several ones, each one corresponding to one single frequency)
ai=context.get_input_parameters(0)
