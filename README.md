# Adaptive-Smearing: This is a Fortran90 code to compute the adaptive-smearing based density of states (DOS).
# After executed the program, you will get the resulting data files: dos_Bloch.dat and dos_Gf.dat, those are the total DOS in the bLoch basis and Wannier basis. 
# You may find orbital-resolved DOS data file as well. 
# The H0_tb.dat file contains the tight-binding parametrs of SrVO3 material obtained from WIEN2k and Wannier90.

You may find ksum.f90 file, which needs to compiled as below:
# mpif90 -o ksum.out ksum.f90 lapack_cmatinv.f

You can run the executable "ksum.out" in one or more number of processors
# mpirun -np 10 ./ksum.out 

You may have asked to provide three input parameters (broadening factor, smearing, NMESH). You can either give the suggested values or different values as well.
