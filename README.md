laplace_fortran
===============

Fortran code using fast integral equation methods to solve Laplace's equation

To play around with FMM only, compile the following files:
 fmm_driver.f
 dapif*.f
 prini.f
 random.f

The directory GMRES contains all of the files needed for GMRES
It would probably be best to compile these into a library

Other files that you may or may not want to take a look at:
 driver_random.f - a driver routine for random.f
