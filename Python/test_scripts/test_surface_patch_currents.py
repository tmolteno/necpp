#instructions to compile :
#
#swig -c++ -python PyNEC.i 
#g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#g++ -shared -lstdc++ nec_context.o nec_output.o c_plot_card.o c_geometry.o misc.o nec_exception.o nec_ground.o c_ggrid.o matrix_algebra.o nec_radiation_pattern.o nec_structure_currents.o c_evlcom.o PyNEC_wrap.o -o _PyNEC.so

#example6.nec
#
#CECYLINDER WITH ATTACHED WIRES
#SP 0 0 10 0 7.3333 0. 0. 38.4
#SP 0 0 10 0 0. 0. 0. 38.4
#SP 0 0 10 0 -7.3333 0. 0. 38.4
#GM 0 1 0. 0. 30.
#SP 0 0 6.89 0. 11. 90. 0. 44.88
#SP 0 0 6.89 0. -11. -90. 0. 44.88
#GR 0 6
#SP 0 0 0. 0. 11. 90. 0. 44.89
#SP 0 0 0. 0. -11. -90. 0. 44.89
#GW 1 4 0. 0. 11. 0. 0. 23. .1 
#GW 2 5 10. 0. 0. 27.6 0. 0. .2 
#GS 0 0 .01
#GE
#FR 0 1 0 0 465.84
#CP 1 1 2 1
#EX 0 1 1 0 1.
#RP 0 73 1 1000 0. 0. 5. 0.
#EX 0 2 1 0 1.
#XQ
#EN

from PyNEC import *

#creation of a nec context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add a patch to the geometry
geo.arbitrary_shaped_patch(10, 0, 7.3333, 0., 0., 38.4)

#add a patch to the geometry
geo.arbitrary_shaped_patch(10, 0, 0, 0., 0., 38.4)

#add a patch to the geometry
geo.arbitrary_shaped_patch(10, 0, -7.3333, 0., 0., 38.4)

#move the structure (here the structure is copied, its copy is rotated by 30 degrees about Z-axis)
geo.move(0, 0, 30, 0, 0, 0, 0, 1, 0)

#add a patch to the geometry
geo.arbitrary_shaped_patch(6.89, 0., 11., 90., 0., 44.88)

#add a patch to the geometry
geo.arbitrary_shaped_patch(6.89, 0., -11., -90., 0., 44.88)

#ask for a cylindrical structure to be generated from the existing structure
geo.generate_cylindrical_structure(0, 6)

#add a patch to the geometry
geo.arbitrary_shaped_patch(0, 0, 11, 90, 0, 44.89)

#add a patch to the geometry
geo.arbitrary_shaped_patch(0, 0, -11, -90, 0, 44.89)

#add a wire to the geometry
geo.wire(1, 4, 0, 0, 11, 0, 0, 23, .1, 1, 1)

#add a wire to the geometry
geo.wire(2, 5, 10, 0, 0, 27.6, 0, 0, .2, 1, 1)

#scale all the structure dimensions by a constant
geo.scale(0.01)

#end of the geometry input
context.geometry_complete(0)

#add a "fr" card to specify the frequency 
context.fr_card(0, 1, 465.84e6, 0)

#add a "ex" card to specify an excitation
context.ex_card(0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(0, 73, 1, 1, 0, 0, 0, 0, 0, 5, 0, 0, 0)

#add a "ex" card to specify an excitation
context.ex_card(0, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0)

#add a "xq" card  to force the simulation execution
context.xq_card(0)

#get the currents
sc = context.get_structure_currents(0)
