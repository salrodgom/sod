## This makefile must be executed with gmake (gnu make).
f90comp = gfortran
FFLAGS ?= -O3 -march=native -funroll-loops -Wall -funroll-loops -fno-protect-parens -ffast-math -ffree-line-length-none -std=f2008
OMPFLAGS ?= -fopenmp

all: sod

sod: src/sod.f90 src/spbe.f90 src/comb.f90 src/canonical.f90 src/compare.f90 src/eqmatrix.f90 src/gc.f90 src/mc.f90 src/exact.f90 src/setup.f90 src/clean.f90 src/entropy.f90 src/base.f90 src/settings.f90 src/energy_calculations.f90 src/core.f90 src/inputs.f90 src/configurations.f90 src/structure_io.f90 src/cell.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/base.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/settings.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/inputs.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/clean.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/configurations.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/symmetry.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/cell.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/structure_io.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/core.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calculations.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/workspace.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/calibration.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/comb.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/setup.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/mc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/exact.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/entropy.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/gc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/canonical.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/compare.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/spbe.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/eqmatrix.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod sod.o eqmatrix.o spbe.o compare.o canonical.o gc.o comb.o entropy.o mc.o exact.o setup.o calibration.o workspace.o energy_calculations.o core.o structure_io.o cell.o symmetry.o configurations.o clean.o inputs.o settings.o base.o
	rm -f *.o *.mod

clean:
	rm -rf bin/sod
