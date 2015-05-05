#!/bin/sh
# cd interface_files
# swig -c++ -python PyNEC.i
# g++ -c ../src/nec_context.cpp PyNEC_wrap.cxx -I../../src -I/usr/local/include/python2.4 -I/usr/local/lib/python2.4/config -DHAVE_CONFIG_H
#
# New script here
#
PYTHON=python2.7
PYDIR=/usr/
NEC_SRC=../src
BUILD=build
rm -rf ${BUILD}
mkdir -p ${BUILD}
cp ${NEC_SRC}/* ${BUILD}
cp ../config.h ${BUILD}
cp interface_files/* ${BUILD}
cd ${BUILD}
swig -c++ -python PyNEC.i
g++ -c nec_context.cpp PyNEC_wrap.cxx -I/usr/lib/python2.7/dist-packages/numpy/numarray/include/ -I/usr/include/python2.7 -I/usr/lib/python2.7/config -DHAVE_CONFIG_H