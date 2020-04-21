#!/bin/bash

gfortran -o premix/dmath.o -c premix/dmath.f -fPIC -fno-range-check
gfortran -o premix/twopnt.o -c premix/twopnt.f
gfortran -o premix/tranlib.o -c premix/tranlib.f
gfortran -o premix/cklib.o -c premix/cklib.f
gfortran -o premix/premix.o -c premix/premix.f
gfortran -c main.f90

gfortran -o get_transport premix/premix.o main.o premix/cklib.o premix/tranlib.o premix/twopnt.o premix/dmath.o