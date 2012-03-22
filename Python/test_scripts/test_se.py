#instructions to compile :
#
#swig -c++ -python PyNEC.i 
#g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#g++ -shared -lstdc++ nec_context.o nec_output.o c_plot_card.o c_geometry.o misc.o nec_exception.o nec_ground.o c_ggrid.o matrix_algebra.o nec_radiation_pattern.o nec_structure_currents.o c_evlcom.o PyNEC_wrap.o -o _PyNEC.so

#example5.nec (modified because there was a little bug in the last TL card...) 
#
#CM 12 ELEMENT LOG PERIODIC ANTENNA IN FREE SPACE 
#CM 78 SEGMENTS. SIGMA=O/L RECEIVING AND TRANS. PATTERNS. 
#CM DIPOLE LENGTH TO DIAMETER RATIO=150. 
#CE TAU=0.93. SIGMA=0.70. BOOM IMPEDANCE=50. OHMS. 
#GW 1 5 0.0000 -1.0000 0.0000000 0.00000 1.0000 0.000 .00667 
#GW 2 5 -.7527 -1.0753 0. -.7527 1.0753 0. .00717 
#GW 3 5 -1.562 -1.1562 0. -1.562 1.1562 0. .00771 
#GW 4 5 -2.4323 -1.2432 0. -2.4323 1.2432 0. .00829 
#GW 5 5 -3.368 -1.3368 0. -3.368 1.3368 0. .00891 
#GW 6 7 -4.3742 -1.4374 0. -4.3742 1.4374 0. .00958 
#GW 7 7 -5.4562 -1.5456 0. -5.4562 1.5456 0. .0103 
#GW 8 7 -6.6195 -1.6619 0. -6.6195 1.6619 0. .01108 
#GW 9 7 -7.8705 -1.787 0. -7.8705 1.787 0. .01191 
#GW 10 7 -9.2156 -1.9215 0. -9.2156 1.9215 0. .01281 
#GW 11 9 -10.6619 -2.0662 0. -10.6619 2.0662 0. .01377 
#GW 12 9 -12.2171 -2.2217 0. -12.2171 2.2217 0. .01481 
#GE 
#FR 0 0 0 0 46.29 0. 
#TL 1 3 2 3 -50. 
#TL 2 3 3 3 -50. 
#TL 3 3 4 3 -50. 
#TL 4 3 5 3 -50. 
#TL 5 3 6 4 -50. 
#TL 6 4 7 4 -50. 
#TL 7 4 8 4 -50. 
#TL 8 4 9 4 -50. 
#TL 9 4 10 4 -50. 
#TL 10 4 11 5 -50. 
#TL 11 5 12 5 -50. 0. 0. 0. .02 
#EX 0 1 3 10 1 0
#RP 0 37 1 1110 90. 0. -5. 0.
#EN

from PyNEC import *

#creation of a nec_context
context=nec_context()

#get the associated geometry
geo = context.get_geometry()

#add wires to the geometry 
geo.wire(1, 5, 0, -1, 0, 0, 1, 0, .00667, 1, 1) 
geo.wire(2, 5, -.7527, -1.0753,0, -.7527, 1.0753, 0, .00717, 1, 1) 
geo.wire(3, 5, -1.562, -1.1562, 0, -1.562, 1.1562, 0, .00771, 1, 1) 
geo.wire(4, 5, -2.4323, -1.2432, 0, -2.4323, 1.2432, 0,.00829, 1, 1) 
geo.wire(5, 5, -3.368, -1.3368, 0, -3.368, 1.3368, 0, .00891, 1, 1) 
geo.wire(6, 7, -4.3742, -1.4374, 0, -4.3742, 1.4374, 0, .00958, 1, 1) 
geo.wire(7, 7, -5.4562, -1.5456, 0, -5.4562, 1.5456, 0, .0103, 1, 1) 
geo.wire(8, 7, -6.6195, -1.6619, 0, -6.6195, 1.6619, 0, .01108, 1, 1) 
geo.wire(9, 7, -7.8705, -1.787, 0, -7.8705, 1.787, 0, .01191, 1, 1) 
geo.wire(10, 7, -9.2156, -1.9215, 0, -9.2156, 1.9215, 0, .01281, 1, 1) 
geo.wire(11, 9, -10.6619, -2.0662, 0, -10.6619, 2.0662, 0, .01377, 1, 1) 
geo.wire(12, 9, -12.2171, -2.2217, 0, -12.2171, 2.2217, 0,.01481, 1, 1)

#end of the geometry input
context.geometry_complete(0)

#add a "fr" card to specify the frequency 
context.fr_card(0, 0, 46.29e6, 0)

#add "tl" cards to generate transmission lines
context.tl_card(1, 3, 2, 3, -50, 0.7527, 0, 0, 0, 0)
context.tl_card(2, 3, 3, 3, -50, 0.8093, 0, 0, 0, 0)
context.tl_card(3, 3, 4, 3, -50, 0.8703, 0, 0, 0, 0)
context.tl_card(4, 3, 5, 3, -50, 0.9357, 0, 0, 0, 0)
context.tl_card(5 ,3 ,6, 4, -50, 1.0062, 0, 0, 0, 0)
context.tl_card(6, 4, 7, 4, -50, 1.082, 0, 0, 0, 0)
context.tl_card(7, 4, 8, 4, -50, 1.1633, 0, 0, 0, 0)
context.tl_card(8, 4, 9, 4, -50, 1.251, 0, 0, 0, 0)
context.tl_card(9, 4, 10, 4, -50, 1.3451, 0, 0, 0, 0)
context.tl_card(10, 4 ,11 , 5, -50, 1.4463, 0, 0, 0, 0)
context.tl_card(11, 5, 12, 5, -50, 1.5552, 0, 0, 0, 0.02)

#add a "ex" card to specify an excitation
context.ex_card(0, 1, 3, 1, 0, 1, 0, 0, 0, 0, 0)

#add a "rp" card to specify radiation pattern sampling parameters and to cause program execution
context.rp_card(0, 37, 1, 1, 1, 1, 0, 90, 0, -5, 0, 0, 0)

#get structure excitation
se = context.get_structure_excitation(0)
