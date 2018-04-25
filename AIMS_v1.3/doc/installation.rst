Installation
============

As of version 1.2, a few strategic parts of the code have been rewritten in
FORTRAN thus leading to a considerable speed up.  These FORTRAN subroutines
are then integrated into the AIMS code thanks to the
`f2py <https://github.com/pearu/f2py/wiki>`_ project.  Accordingly, these
FORTRAN subroutines need to be compiled before running AIMS.  A Makefile
has been provided for convenience.  Hence, one simply needs to type the
command::

    make

The user may change the choice of FORTRAN compiler as well as the compilation
options by editing the Makefile.
