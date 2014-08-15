import sys
import os

print('Installation start')
#==================================================================================================================================================================================
#get the platform used
platform = sys.platform

#==================================================================================================================================================================================
#checks the Numarray installation

print('\nNumarray installation check...')
try:
	from numpy import numarray
except:
	print('FAILURE - PyNEC requires Numarray to be installed. Please install Numarrray first.')
	print('Installation aborted')
	sys.exit(-1)

print('SUCCESS - Numarray seems to be installed.')

#==================================================================================================================================================================================
#Detects the Python version number

print('\nDetecting Python version number...')
print(sys.version)
python_version = sys.version[:3]
if python_version != '2.3' and python_version != '2.4' :
	print('WARNING - This python version has not yet been tested. Try and move to Python version 2.3 or 2.4 if you encounter problems')
api_version = sys.api_version
if api_version != 1012 :
	print('WARNING - This api version has not yet been tested. Try and move to API version 1012 if you encounter problems')

#==================================================================================================================================================================================
#Looks for the python "include path"

print('\nLooking for the Python include path...')

def test_include_path(str):
	open(os.path.join(str, 'Python.h'), 'rb')
	return 0

def search_include_path(flag):
	try:
		if platform == 'win32' :
			include_path = os.path.join(sys.prefix, 'include')
		else :
			include_path = os.path.join(sys.prefix, 'include', 'python'+python_version)

		test_include_path(include_path)
		return include_path
	except:
		if flag == False :
                       print("PROBLEM - Can not find the Python include path.\nPlease enter the actual Python include path ( example : " + include_path + ") :")
		else :
                       print('Try again :')

		include_provided = False
		while(include_provided == False) :
			include_path = sys.stdin.readline()[:-1]
			include_provided = True
		try :
			test_include_path(include_path)
			return include_path
		except :
			print('Wrong Python include path')
			return search_include_path(True)

		if flag == False :
			print('The include path provided seems to be correct')


python_include_path = search_include_path(False)
print('SUCCESS - Python include path found')

#==================================================================================================================================================================================
#Looks for the Numarray "API path"

print('\nLooking for the Numarray API path...')

def test_api_path(str):
	open(os.path.join(str, 'libnumarray.h'), 'rb')
	return 0

def search_api_path(flag):
	try:
		if platform == 'win32' :
			api_path = os.path.join(sys.prefix, 'include','numarray')
		else :
			api_path = os.path.join(sys.prefix, 'include', 'python'+python_version, 'numarray')
		test_api_path(api_path)
		return api_path
	except:
		if flag == False :
			print("PROBLEM - Can not find the Numarray API path.\nPlease enter the actual Numarray API path ( example : " + api_path + ") :")
		else :
			print('Try again :')

		api_provided = False
		while(api_provided == False) :
			api_path = sys.stdin.readline()[:-1]
			api_provided = True
		try :
			test_api_path(api_path)
			return api_path
		except :
			print('Wrong Numarray API path')
			return search_api_path(True)

		if flag == False :
			print('The API path provided seems to be correct')


numarray_api_path = search_api_path(False)
print('SUCCESS - Numarray API path found')

#==================================================================================================================================================================================
#Looks for the python "lib path"

print('\nLooking for the Python lib path...')

def test_lib_path(str):
	open(os.path.join(str, 'os.py'), 'rb')
	return 0

def search_lib_path(flag):
 	try:
		if platform == 'win32' :
			lib_path = os.path.join(sys.prefix, 'lib')
		else :
			lib_path = os.path.join(sys.prefix, 'lib', 'python'+python_version)
		test_lib_path(lib_path)
		return lib_path
	except:
		if flag == False :
			print("PROBLEM - Can not find the Python lib path.\nPlease enter the actual Python lib path ( example : " + lib_path + ") :")
		else :
			print('Try again :')

		lib_provided = False
		while(lib_provided == False) :
			lib_path = sys.stdin.readline()[:-1]
			lib_provided = True
		try :
			test_lib_path(lib_path)
			return lib_path
		except :
			print('Wrong Python lib path')
			return search_lib_path(True)

		if flag == False :
			print('The lib path provided seems to be correct.')


python_lib_path = search_lib_path(False)
print('SUCCESS - Python lib path found')

#==================================================================================================================================================================================
#Generates the makefile

print('\nGenerating the makefile...')

if platform == 'linux2' or platform == 'win32' :
	try :
		f = open ('makefile.'+platform, 'wb')
	except :
		print('FAILURE - the install process failed to create a file.')
		print('Installation aborted')
		sys.exit(-1)

	print('File created')

       	print('Writting the file...')
       	try :
		f.write('# MAKEFILE FOR THE PYTHON WRAPPING OF NEC-2\n')
		
		l = ['misc', 'c_ggrid', 'c_evlcom', 'c_plot_card', 'matrix_algebra', 'nec_output', 'c_geometry']
		list_hdrs = l + ['math_util', 'common', 'nec_results']
		list_src = l + ['nec_context', 'nec_ground', 'nec_radiation_pattern', 'nec_exception', 'nec_structure_currents']
		
		if platform == 'linux2' :
			f.write('\nPYMOD=python_module/_PyNEC.so\n');
			
			hdrs = '\nHDRS 	= '
			i = -1
			for j in range(list_hdrs.__len__()) :
				i += 1
				if i == 4 :
					hdrs += "\\\n	"
					i = 0
				hdrs += "wrap/src/" + list_hdrs[j] + ".h" + " " 
			f.write(hdrs + "\n")
			
			src = '\nSRC 	= '
			i = -1
			for j in range(list_src.__len__()) :
				i += 1
				if i == 4 :
					src += "\\\n	"
					i = 0
				src += "wrap/src/" + list_src[j] + ".cpp" + " "
			f.write(src + "\n")
			
			f.write("\nOBJS 	= $(SRC:.cpp=.o)\n")

			f.write("\nCXX=g++\n")
			
			f.write("\nLDFLAGS=-shared -lstdc++\n") 

			f.write("\nCXXFLAGS=-Wall\n")
			
			swigcxxflags = "\nSWIGCXXFLAGS=-c"
			swigcxxflags += " -I"+python_include_path
			swigcxxflags += " -I"+numarray_api_path
			swigcxxflags += " -I"+python_lib_path
			swigcxxflags += " -DHAVE_CONFIG_H"
			f.write(swigcxxflags+"\n")
			
			f.write("\nSWIGFLAGS=-python -c++\n")
			
			f.write("\n\n")
			
			f.write("\nall: $(PYMOD)\n")
			
			pynec_so = "\npython_module/_PyNEC.so: $(OBJS) wrap/PyNEC_wrap.o\n"
			pynec_so += "	cd ..\n"
			pynec_so += "	$(CXX) $(LDFLAGS) $^ -o $@\n" 
			pynec_so += "	@echo\n"
			pynec_so += "	@echo To install the module on your python system, type make -f makefile." + platform + " install - root privileges may be required"
			f.write(pynec_so)
			
			pynec_wrap_o = "\nwrap/PyNEC_wrap.o: wrap/PyNEC_wrap.cxx\n"
			pynec_wrap_o += "	$(CXX) $(SWIGCXXFLAGS) $< wrap/src/nec_context.cpp\n"
			pynec_wrap_o += "	@mv PyNEC_wrap.o wrap/\n"
			pynec_wrap_o += "	@rm -rf nec_context.o\n"
			f.write(pynec_wrap_o)
			
			pynec_wrap_cxx = "\nwrap/PyNEC_wrap.cxx: wrap/PyNEC.i\n"
			pynec_wrap_cxx += "	swig $(SWIGFLAGS) $<\n"
			pynec_wrap_cxx += "	@rm -rf wrap/PyNEC.py\n"
			f.write(pynec_wrap_cxx)
			
			o_from_cpp = "\n%.o: %.cpp\n"
			o_from_cpp += "	$(CXX) $(CXXFLAGS) -o $@ -c $<\n"
			f.write(o_from_cpp)
			
			install = "\ninstall:\n"
			install += "	@rm -fr " + python_lib_path + "/site-packages/PyNEC\n"
			install += "	@mkdir " + python_lib_path + "/site-packages/PyNEC\n"
			install += "	cp python_module/*.* " + python_lib_path + "/site-packages/PyNEC/\n"
			f.write(install)
			
			clean = "\nclean:\n"
			clean += "	rm -rf wrap/src/*.o wrap/PyNEC_wrap.cxx wrap/PyNEC_wrap.o python_module/*.pyc python_module/_PyNEC.so"
			f.write(clean)
			
			print('File written\n')
			print('You can now use the generated makefile (named makefile.' + platform +') to compile PyNEC\n')
	       
	except :
		print('FAILURE - an error occured so the file has not been filled properly.')
		print('Installation aborted')
		sys.exit(-1)
else :
	print('FAILURE - The installation of PyNEC has not yet been tested on this platform.')
	print('Installation aborted')
	sys.exit(-1)

