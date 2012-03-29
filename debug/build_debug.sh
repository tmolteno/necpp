#!/bin/sh
../configure --prefix=`pwd` CPPFLAGS="-DDEBUG -DNEC_ERROR_CHECK" CXXFLAGS="-g -O0"
make
make install

