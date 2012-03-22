#!/bin/sh
gcc -Wall -o atlas_check atlas_check.cpp -L /usr/lib/atlas -llapack -lblas -latlas -lstdc++
./atlas_check
rm atlas_check
gfortran -Wall -o atlas_check atlas_check.f
./atlas_check
rm atlas_check
gfortran -Wall -o atlas_check lu_ralston.f
./atlas_check
rm atlas_check
