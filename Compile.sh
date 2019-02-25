#! /bin/bash -x
rm *.o *.mod *.exe
gfortran -c Constants.f90 
gfortran -c BARSumRoutines.f90
gfortran -o FEP-OPT.exe BARsum.f90 BARSumRoutines.o Constants.o
rm *.o *.mod 
