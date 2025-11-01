## This makefile must be executed with gmake (gnu make).

f90comp = gfortran
OMPFLAGS ?= -fopenmp

all:
	$(f90comp) -o bin/combsod src/factorials.f90 src/bubble.f90  src/ksubset.f90  src/member.f90  src/cell.f90 src/ccf.f90 src/combsod.f90
	$(f90comp) -o bin/genersod src/member.f90 src/cell.f90 src/genersod.f90
	$(f90comp) -o bin/spbesod src/spbesod.f90
	$(f90comp) -o bin/invertOUTSOD src/invertOUTSOD.f90
	$(f90comp) -o bin/statsod  src/statsod.f90
	$(f90comp) -o bin/gcstatsod  src/factorials.f90 src/momenta.f90 src/gcstatsod.f90
	$(f90comp) -o bin/peaks2spec  src/peaks2spec.f90
	$(f90comp) $(OMPFLAGS) -c src/energy_calc.f90
	$(f90comp) $(OMPFLAGS) -o bin/sod_boltzmann_mc src/sod_boltzmann_mc.f90 energy_calc.o
#	$(f90comp) -c src/mc_sampler.f90
#	$(f90comp) -o bin/sod_mc energy_calc.o mc_sampler.o src/sod_mc.f90
	rm -f *.o *.mod

sod_boltzmann_mc: src/sod_boltzmann_mc.f90 src/energy_calc.f90
	$(f90comp) $(OMPFLAGS) -c src/energy_calc.f90
	$(f90comp) $(OMPFLAGS) -o bin/sod_boltzmann_mc src/sod_boltzmann_mc.f90 energy_calc.o
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
	rm bin/sod_mc


