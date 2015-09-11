#!/bin/bash
#
#  Compiles fracd.f90 with FFTW
#
F90 = gfortran
FFTWFLAGS = -L/usr/local/lib -I/usr/local/include -lfftw3

fracd: fracd.f90
	$(F90) -o fracd fracd.f90 $(FFTWFLAGS)
