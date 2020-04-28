#!/bin/bash

gfortran -g -o src/chemkin_m.o -c src/chemkin_m.f90 -J src
gfortran -g -o src/dmath.o     -c src/dmath.f -fPIC -fno-range-check
gfortran -g -o src/twopnt.o    -c src/twopnt.f
gfortran -g -o src/tranlib.o   -c src/tranlib.f
gfortran -g -o src/cklib.o     -c src/cklib.f
gfortran -g -o src/dasac.o     -c src/dasac.f 
gfortran -g -o src/senkin.o    -c src/senkin.f -I src
gfortran -g -o src/main.o      -c src/main.f90 -I src

gfortran -g -o get_CFD_values src/main.o src/chemkin_m.o src/dasac.f src/senkin.o src/cklib.o src/tranlib.o src/dmath.o