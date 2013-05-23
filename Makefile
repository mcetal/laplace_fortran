#laplace_fortran/Makefile

OBJECTS = fmm_driver.f dapif*.f prini.f random.f

.PHONY: test

fmmtest.exe: $(OBJECTS)
	gfortran -fno-automatic -std=legacy $(OBJECTS) -o fmmtest.exe

test: fmmtest.exe
	@echo Testing FMM...
	./fmmtest.exe


