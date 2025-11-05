## This makefile must be executed with gmake (gnu make).

f90comp = gfortran
FFLAGS ?= -O3 -march=native -funroll-loops
OMPFLAGS ?= -fopenmp

all:
	$(f90comp) $(FFLAGS) -o bin/combsod src/factorials.f90 src/bubble.f90 src/ksubset.f90 src/member.f90 src/cell.f90 src/ccf.f90 src/combsod.f90
	$(f90comp) $(FFLAGS) -o bin/genersod src/member.f90 src/cell.f90 src/genersod.f90
	$(f90comp) $(FFLAGS) -o bin/spbesod src/spbesod.f90
	$(f90comp) $(FFLAGS) -o bin/invertOUTSOD src/invertOUTSOD.f90
	$(f90comp) $(FFLAGS) -o bin/statsod src/statsod.f90
	$(f90comp) $(FFLAGS) -o bin/gcstatsod src/factorials.f90 src/momenta.f90 src/gcstatsod.f90
	$(f90comp) $(FFLAGS) -o bin/peaks2spec src/peaks2spec.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_boltzmann_base.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_boltzmann_mc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_boltzmann_mc sod_boltzmann_mc.o sod_boltzmann_base.o energy_calc.o
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_boltzmann_exact src/sod_boltzmann_exact.f90 sod_boltzmann_base.o energy_calc.o
#	$(f90comp) -c src/mc_sampler.f90
#	$(f90comp) -o bin/sod_mc energy_calc.o mc_sampler.o src/sod_mc.f90
	rm -f *.o *.mod

sod_boltzmann_mc: src/sod_boltzmann_mc.f90 src/sod_boltzmann_base.f90 src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_boltzmann_base.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_boltzmann_mc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_boltzmann_mc sod_boltzmann_mc.o sod_boltzmann_base.o energy_calc.o
	rm -f *.o *.mod

sod_boltzmann_exact: src/sod_boltzmann_exact.f90 src/sod_boltzmann_mc.f90 src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_boltzmann_base.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_boltzmann_exact src/sod_boltzmann_exact.f90 sod_boltzmann_base.o energy_calc.o
	rm -f *.o *.mod

clean:
	rm bin/combsod  
	rm bin/invertOUTSOD  
	rm bin/statsod
	rm bin/gcstatsod
	rm bin/genersod
	rm bin/spbesod
	rm bin/peaks2spec
	rm bin/sod_boltzmann_mc
	rm bin/sod_boltzmann_exact


