F90 = gfortran
FLG = -O2 #-fcheck=all #-fbounds-check

aims_fortran: aims_fortran.f90
	f2py -c aims_fortran.f90 -m aims_fortran --fcompiler=$(F90) --f90flags=$(FLG)

clean:
	rm -f aims_fortran.so

