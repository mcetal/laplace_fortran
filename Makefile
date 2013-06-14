#laplace_fortran/Makefile

OBJECTS1 = fmm_driver.f dapif*.f prini.f random.f
OBJECTS2 = dapif2.f dcfft.f matplot.f random.f targets.m dapif1.f dapif3.f laplace_driver.f prini.f   
LIB =   libhfmm2d.a GMRES/*.f
FC  = gfortran
FLAGS = 

.PHONY: test clean

fmmtest.exe: $(OBJECTS)
	$(FC) -fno-automatic -std=legacy $(OBJECTS) -o fmmtest.exe

test: fmmtest.exe
	@echo Testing FMM...
	./fmmtest.exe
new: new.exe
	@echo Checking stuff
	./new.exe
new.exe: new.f
	$(FC) -fno-automatic -std=legacy new.f -o new.exe
clean: 
	rm *.exe *.o

