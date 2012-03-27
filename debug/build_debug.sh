#!/bin/sh
../configure --prefix=/dbg    CPPFLAGS="-DDEBUG -DNEC_ERROR_CHECK" CXXFLAGS="-g -O0" && make && make install

