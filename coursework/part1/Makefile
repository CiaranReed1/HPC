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
SPEEDFILE ?= speedtest_1.datS
SRC= ?= part1.c
.PHONY: scaling
scaling:
	@echo "Submitting strong scaling jobs..."
	@echo "cores,time_seconds" > $(SPEEDFILE)
	@for c in $(CORES); do \
		echo "  -> Submitting job with $$c cores"; \
		sbatch -c $$c --export=ALL,SPEEDFILE=$(SPEEDFILE),SRC=$(SRC) speedtest1.sbatch; \
	done
