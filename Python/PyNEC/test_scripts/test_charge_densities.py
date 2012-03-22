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

print('beginning of the test')
import sys

print('import of the module')
from PyNEC import *

print('beginning of card input')
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

print('end of card input\n')
print('get the currents (here the currents and charge densities datas)\n')
sc = context.get_structure_currents(0)

print('get the array of length of the segment wires for the printing of currents')
print(sc.get_current_segment_length())

print('\nget the array of y-coordinates of segment centers for the printing of charge densities')
print(sc.get_q_density_segment_center_y())

print('\ntry and get the array of x-coordinates of patch centers')
try:
	print(sc.get_patch_center_x())
	print('ERROR - there are no patches !')
	sys.exit(-1)
except Warning, msg :
	print('Warning : '+msg.__str__())
	print("SUCCESS - an exception is raised as there are no patches their datas can't be requested")

print('\nend of the test - you can compare the results with the ones provided by NEC-2 using test_charge_densities.nec as the input file\n')
