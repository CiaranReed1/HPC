# Simple cleanup Makefile

# List of binaries to remove (add names if needed)
BINS = serial part1 

# SLURM output files
SLURM_OUT = slurm-*.out *.err

# Other generated files
DATA = serial.dat part1.dat correctness_check.txt

.PHONY: clean

clean:
	@echo "Cleaning binaries, SLURM outputs, and data files..."
	@rm -f $(BINS) $(SLURM_OUT) $(DATA)

CORES = 1 2 4 8 16 32 64 128

.PHONY: scaling
scaling:
	@echo "Submitting strong scaling jobs..."
	@echo "cores,time_seconds" > speedtest_1.dat
	@for c in $(CORES); do \
		echo "  -> Submitting job with $$c cores"; \
		sbatch -c $$c speedtest1.sbatch; \
	done