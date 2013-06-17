laplace_fortran : A Laplace Equation solver written in fortran
==============================================================

Fortran code using fast integral equation methods to solve Laplace's equation

To play around with the [Fast Multipole Method](http://en.wikipedia.org/wiki/Fast_multipole_method) : 
    
    make test

To run the laplace equation solver:
    
    make laplace

To test any small tweaks, put them in new.f and use
    
    make new

Use
    make help
when needed.

For FMM, The following files are used:
* fmm_driver.f
* dapif*.f(old fmm files)
* prini.f
* random.f

    
    libhfmm2d.a 

is a compiled library of a [new fmm code](http://www.cims.nyu.edu/cmcl/fmm2dlib/fmmlib2d-1.2.zip) for Laplace and Helmholtz equation.
To get the library,
    
    #Extract the files
    make lib

The directory GMRES contains all of the files needed for GMRES
It would probably be best to compile these into a library

Other files that you may or may not want to take a look at:
 driver_random.f - a driver routine for random.f


