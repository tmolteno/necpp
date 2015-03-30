# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
setup.py file for necpp Python module
"""

from distutils.core import setup, Extension
from glob import glob
import os

nec_sources = ['necpp.i']
nec_sources = [fn for fn in glob('../src/*.cpp') 
         if not os.path.basename(fn).endswith('_tb.cpp')
         if not os.path.basename(fn).startswith('nec2cpp.cpp')
         if not os.path.basename(fn).startswith('necDiff.cpp')]
nec_sources.extend(glob("necpp_wrap.c"))

nec_headers = []
nec_headers.extend(glob("../src/*.h"))
nec_headers.extend(glob("../config.h"))


# At the moment, the config.h file is needed, and this should be generated from the ./configure
# command in the parent directory. Use ./configure --without-lapack to avoid dependance on LAPACK
#
necpp_module = Extension('_necpp',
    sources=nec_sources,
    include_dirs=['../src/', '../'],
    depends=nec_headers,
    define_macros=[('BUILD_PYTHON', '1')],
    language='c++'
    )

# 
#   
#necpp_module = Extension('_necpp',
                           #sources=['necpp_wrap.c'],
                           #include_dirs=['/usr/local/include'],
                           #libraries=['necpp']
                           #)


setup (name = 'necpp',
       version = '0.1.4',
       author  = "Tim Molteno",
       author_email  = "tim@physics.otago.ac.nz",
       url  = "http://github.com/tmolteno/necpp",
       description = "Python Antenna Simulation Module (nec2++)",
       ext_modules = [necpp_module],
       dependency_links=['http://github.com/tmolteno/necpp'],
       py_modules = ["necpp"],
       license='GPLv2'
       )
