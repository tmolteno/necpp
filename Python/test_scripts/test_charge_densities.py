#instructions to compile :
#
#swig -c++ -python PyNEC.i 
#g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#g++ -shared -lstdc++ nec_context.o nec_output.o c_plot_card.o c_geometry.o misc.o nec_exception.o nec_ground.o c_ggrid.o matrix_algebra.o nec_radiation_pattern.o nec_structure_currents.o c_evlcom.o PyNEC_wrap.o -o _PyNEC.so

#dipole_anim.nec
#
#CM Simple dipole, with calculation of currents, charges and near field.
#CE
#GW  1   21         0      -0.25      0.0         0       0.25       0.0     0.001
#GE
#EX  0   1   11   00         1         0
#PQ  0,  0
#NE  0,  1,20,20,  0,0.05,0.05,  0,0.05,0.05
#NH  0,  1,20,20,  0,0.05,0.05,  0,0.05,0.05
#RP  0, 19, 36, 1000, 0, 0, 10, 10
#EN

from PyNEC import *

#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add a wire to the geometry
geo.wire(1, 21, 0, -0.25, 0.0, 0, 0.25, 0.0, 0.001, 1, 1)

#end of the geometry input
context.geometry_complete(0)

#add a "ex" card to specify an excitation
context.ex_card(0, 1, 11, 0, 0, 1, 0, 0, 0, 0, 0)

#add a "pq" card to ask for the charge densities to be computed
context.pq_card(0, 0, 0, 0)

#add a "ne" card to ask for the "near electric field pattern" to be computed
context.ne_card(0, 1, 20, 20, 0, 0.05, 0.05, 0, 0.05, 0.05)

#add a "nh" card to ask for the "near magnetic field pattern" to be computed
context.nh_card(0, 1, 20, 20, 0, 0.05, 0.05, 0, 0.05, 0.05)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(0, 19, 36, 1, 0, 0, 0, 0, 0, 10, 10, 0, 0)

#get the currents
sc = context.get_structure_currents(0)
