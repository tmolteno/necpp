#!/usr/bin/env python

"""
setup.py file for necpp Python module
"""

from distutils.core import setup, Extension


necpp_module = Extension('_necpp',
                           sources=['necpp_wrap.c'],
                           include_dirs=['/usr/local/include'],
                           libraries=['necpp']
                           )

setup (name = 'necpp',
       version = '0.1.2',
       author      = "Tim Molteno tim@physics.otago.ac.nz",
       description = """Python Antenna Simulation Module (nec2++)""",
       ext_modules = [necpp_module],
       py_modules = ["necpp"],
       )
