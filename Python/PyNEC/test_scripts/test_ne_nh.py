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

#get the currents
sc = context.get_structure_currents(0)

print('end of card input\n')
print('get the near electric and magnetic field\n')
ne = context.get_near_field_pattern(0)
nh = context.get_near_field_pattern(1)

print('try and get the array of z-coordinates of electric field from the near magnetic field')
try:
	print(nh.get_e_x())
	print('ERROR - the result is a magnetic field !')
	sys.exit(-1)
except Warning, msg :
	print('Warning : '+msg.__str__())
	print('SUCCESS - an exception has been raised as the result is a magnetic field\n')

print('get the array of x-coordinates of magnetic field from the near magnetic field')
print(nh.get_h_x())

print('\nend of the test - you can compare the results with the ones provided by NEC-2 using test_ne_nh.nec as the input file\n')
