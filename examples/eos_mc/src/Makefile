FC = gfortran
FFLAGS ?= -O2 -Wall -Wextra -std=f2008
BINDIR ?= ..

TARGETS := $(BINDIR)/sod_boltzmann_mc $(BINDIR)/sod_mc
OBJECTS_COMMON := energy_calc.o
OBJECTS_MC := mc_sampler.o sod_mc.o
OBJECTS_BOLTZ := sod_boltzmann_mc.o

.PHONY: all clean veryclean

all: $(TARGETS)

$(BINDIR):
	@mkdir -p $(BINDIR)

$(BINDIR)/sod_boltzmann_mc: $(OBJECTS_COMMON) $(OBJECTS_BOLTZ) | $(BINDIR)
	$(FC) $(FFLAGS) $^ -o $@

$(BINDIR)/sod_mc: $(OBJECTS_COMMON) $(OBJECTS_MC) | $(BINDIR)
	$(FC) $(FFLAGS) $^ -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f $(TARGETS)
