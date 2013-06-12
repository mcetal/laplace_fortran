#laplace_fortran/Makefile

OBJECTS = fmm_driver.f dapif*.f prini.f random.f

.PHONY: test clean

fmmtest.exe: $(OBJECTS)
	gfortran -fno-automatic -std=legacy $(OBJECTS) -o fmmtest.exe

test: fmmtest.exe
	@echo Testing FMM...
	./fmmtest.exe
new: new.exe
	@echo Checking stuff
	./new.exe
new.exe: new.f
	gfortran -fno-automatic -std=legacy new.f -o new.exe
clean: 
	rm *.exe *.o

