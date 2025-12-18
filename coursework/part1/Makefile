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