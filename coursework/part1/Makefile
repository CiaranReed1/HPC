# Simple cleanup Makefile

# List of binaries to remove (add names if needed)
BINS = serial part1 

# SLURM output files
SLURM_OUT = *.out *.err

# Other generated files
DATA = part1_*_cores.dat correctness_check_*.txt

.PHONY: clean

clean:
	@echo "Cleaning binaries, SLURM outputs, and data files..."
	@rm -f $(BINS) $(SLURM_OUT) $(DATA)
	@rm bin/*
	@rm outputs/*

CORES = 1 2 4 8 16 32 64 128
SRC= ?= part1.c
BASENAME := $(notdir $(SRC:.c=))
SPEEDFILE := newData/$(BASENAME)_strong_scaling.dat
.PHONY: scaling
scaling:
	@echo "Submitting strong scaling jobs..."
	@echo "cores,time_seconds" > $(SPEEDFILE)
	@for c in $(CORES); do \
		echo "  -> Submitting job with $$c cores"; \
		sbatch -c $$c --export=ALL,SRC=$(SRC),SPEEDFILE=$(SPEEDFILE) speedtest1.sbatch; \
	done
