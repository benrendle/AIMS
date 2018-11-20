F90 = gfortran
FLG = -O2 #-fcheck=all #-fbounds-check
FPY = f2py
#FPY = f2py3 # for python3
#FPY = python3 /full/address/to/f2py # another way of working with python3

aims_fortran: aims_fortran.f90
	#f2py -c aims_fortran.f90 -m aims_fortran --fcompiler=$(F90) --f90flags=$(FLG)
	f2py -c aims_fortran.f90 -m aims_fortran --f90flags=$(FLG)

clean:
	rm -f aims_fortran.so

