#laplace_fortran/Makefile

#source files for fmm testing
FSOURCES = fmm_driver.f dapif1.f dapif2.f dapif3.f prini.f random.f

#source files for solving the laplace equation
LSOURCES = dapif2.f dcfft.f matplot.f random.f dapif1.f dapif3.f laplace_driver.f prini.f 

#object files for fmm testing  
FOBJECTS = $(FSOURCES:.f=.o)

#object files for laplace equation
LOBJECTS = $(LSOURCES:.f=.o) 

#libraries required for laplace equation driver. 
#libfmm2d.a is the compiled library of the latest fmm code.
#GMRES has the required files for the GMRES method 
LIB =   libhfmm2d.a GMRES/*.f

#Set your fortran compiler here
FC  = gfortran

#Flags for fmm testing 
FFLAGS = -fno-automatic -std=legacy
#Flags for laplace equation solver
LFLAGS1 = -c
LFLAGS2 = -o

#type `make clean` to rm the object files
.PHONY: test clean help

fmmtest.exe: $(FOBJECTS)
	$(FC) -fno-automatic -std=legacy $(FSOURCES) -o fmmtest.exe

test: fmmtest.exe
	@echo Testing FMM...
	./fmmtest.exe
new: new.exe
	@echo Checking stuff
	./new.exe
new.exe: new.f
	$(FC) $(FFLAGS) new.f -o new.exe

laplace_driver.exe: $(LSOURCES)
	$(FC) $(LFLAGS1) $(LSOURCES) 
	$(FC) $(LFLAGS2) $@ $(LOBJECTS) $(LIB)  	

laplace: laplace_driver.exe
	@echo Testing Laplace Equation...
	./laplace_driver.exe
clean: 
	rm *.exe $(LOBJECTS) fmm_driver.o

help:
	@echo make test    --- Test FMM
	@echo make new     --- Put your code in new.f and test using this command
	@echo make laplace --- Test Laplace Equation Solver
	@echo make clean   --- Remove all object files 
	@echo make help    --- A self-reference

