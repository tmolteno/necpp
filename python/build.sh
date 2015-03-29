#!/bin/sh
# Script to build the nec2++ ruby module.
# Change the RUBY environment variable 
# as appropriate for your system
../configure --without-lapack
PYTHON=python
swig -v -I/usr/local/include -python necpp.i
sudo python setup.py install
