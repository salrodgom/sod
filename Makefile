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
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_base.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_config.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_symmetry.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_workspace.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_calibration.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_mc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_exact.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_mc_main.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_exact_main.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_ensemble_mc sod_ensemble_mc_main.o sod_ensemble_mc.o sod_ensemble_calibration.o sod_ensemble_workspace.o sod_ensemble_symmetry.o sod_ensemble_config.o sod_ensemble_base.o energy_calc.o
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_ensemble_exact sod_ensemble_exact_main.o sod_ensemble_exact.o sod_ensemble_calibration.o sod_ensemble_symmetry.o sod_ensemble_config.o sod_ensemble_base.o energy_calc.o
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_ensemble sod_ensemble.o sod_ensemble_mc.o sod_ensemble_exact.o sod_ensemble_calibration.o sod_ensemble_workspace.o sod_ensemble_symmetry.o sod_ensemble_config.o sod_ensemble_base.o energy_calc.o
	rm -f *.o *.mod

#sod_ensemble_mc: src/sod_ensemble_mc.f90 src/sod_ensemble_base.f90 src/energy_calc.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calc.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_base.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_config.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_symmetry.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_workspace.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_calibration.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_mc.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_mc_main.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_ensemble_mc sod_ensemble_mc_main.o sod_ensemble_mc.o sod_ensemble_calibration.o sod_ensemble_workspace.o sod_ensemble_symmetry.o sod_ensemble_config.o sod_ensemble_base.o energy_calc.o
#	rm -f *.o *.mod
#
#sod_ensemble_exact: src/sod_ensemble_exact.f90 src/sod_ensemble_base.f90 src/energy_calc.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calc.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_base.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_config.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_symmetry.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_workspace.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_calibration.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_exact.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_exact_main.f90
#	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_ensemble_exact sod_ensemble_exact_main.o sod_ensemble_exact.o sod_ensemble_calibration.o sod_ensemble_symmetry.o sod_ensemble_config.o sod_ensemble_base.o energy_calc.o
#	$(f90comp) $(FFLAGS) -o bin/sod_config_entropy src/sod_ensemble_base.f90 src/sod_config_entropy.f90
#	rm -f *.o *.mod

sod_ensemble: src/sod_ensemble.f90 src/sod_ensemble_mc.f90 src/sod_ensemble_exact.f90 src/sod_ensemble_base.f90 src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/energy_calc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_base.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_config.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_symmetry.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_workspace.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_calibration.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_mc.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble_exact.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -c src/sod_ensemble.f90
	$(f90comp) $(FFLAGS) $(OMPFLAGS) -o bin/sod_ensemble sod_ensemble.o sod_ensemble_mc.o sod_ensemble_exact.o sod_ensemble_calibration.o sod_ensemble_workspace.o sod_ensemble_symmetry.o sod_ensemble_config.o sod_ensemble_base.o energy_calc.o
	rm -f *.o *.mod

sod_config_entropy: src/sod_config_entropy.f90 src/sod_ensemble_base.f90
	$(f90comp) $(FFLAGS) -o bin/sod_config_entropy src/sod_ensemble_base.f90 src/sod_config_entropy.f90
	rm -f *.o *.mod

clean:
	rm bin/combsod  
	rm bin/invertOUTSOD  
	rm bin/statsod
	rm bin/gcstatsod
	rm bin/genersod
	rm bin/spbesod
	rm bin/peaks2spec
	rm bin/sod_ensemble_mc
	rm bin/sod_ensemble_exact
	rm bin/sod_ensemble
	rm bin/sod_config_entropy
