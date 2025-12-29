# =========================
# User-configurable options
# =========================

# Source file to compile (override on command line)
#   make compile MAIN=test.c
MAIN ?= part2.c
SERIAL_SRC ?= serial.c
# MPI runs per job
RUNS ?= 1

# MPI ranks for weak scaling
RANKS = 1 2 4 


# =========================
# Project layout
# =========================

SRC_DIR := src
BIN_DIR := bin
OUT_DIR := outputs
DATA_DIR := scalingData
PROG_DATA_DIR := programData


# =========================
# Derived names
# =========================

BASENAME := $(notdir $(MAIN:.c=))
BIN := $(BIN_DIR)/$(BASENAME)
SERIAL_BIN = bin/serial 

SPEEDFILE := $(DATA_DIR)/$(BASENAME)_weak_scaling.dat


# =========================
# Compiler settings
# =========================

CC := mpicc
CFLAGS := -I.
LDFLAGS := -lm


# =========================
# Phony targets
# =========================

.PHONY: compile scaling clean serial serial-compile


# =========================
# Build rule (single-file)
# =========================

$(BIN): $(SRC_DIR)/$(MAIN) params.h
	@mkdir -p $(BIN_DIR)
	module purge && \
	module load gcc && \
	module load openmpi && \
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@
	@echo "Built $@"

$(SERIAL_BIN): $(SRC_DIR)/$(SERIAL_SRC) params.h
	@mkdir -p $(BIN_DIR)
	gcc -I. $< -lm -o $@
	@echo "Built serial binary: $@"


compile: $(BIN)
serial-compile: $(SERIAL_BIN)

# =========================
# Weak scaling submission
# =========================

scaling: compile
	@mkdir -p $(OUT_DIR) $(DATA_DIR)
	@echo "problem_size,size,run_id,init_time,step_time,dxdt_time,norm_time,total_time" > $(SPEEDFILE)
	@echo "Submitting weak scaling jobs for $(BIN)"
	@for n in $(RANKS); do \
		echo "  -> MPI ranks: $$n"; \
		sbatch -n $$n \
			--export=ALL,PROBLEM_SIZE=$$n,SIZE=$$n,BIN=$(BIN),RUNS=$(RUNS),SPEEDFILE=$(SPEEDFILE) \
			scaling.sbatch; \
		sbatch -n 1 \
			--export=ALL,PROBLEM_SIZE=$$n,SIZE=1,BIN=$(BIN),RUNS=$(RUNS),SPEEDFILE=$(SPEEDFILE) \
			scaling.sbatch; \
	done

#========================
# Serial Submission
#========================
serial: serial-compile
	@mkdir -p $(OUT_DIR) $(DATA_DIR)
	@echo "Submitting serial jobs for $(SERIAL_BIN)"
	@for n in $(RANKS); do \
		echo "  -> Problem Size: $$n"; \
		sbatch \
			--export=ALL,PROBLEM_SIZE=$$n,BIN=$(SERIAL_BIN),RUNS=$(RUNS),SPEEDFILE=$(SPEEDFILE) \
			serial.sbatch; \
	done

# =========================
# Cleanup
# =========================

clean:
	@echo "Cleaning build artifacts and outputs..."
	rm -f $(BIN_DIR)/* \
	      $(OUT_DIR)/*.out $(OUT_DIR)/*.err \
	      $(DATA_DIR)/*.lock
	rm -f $(PROG_DATA_DIR)/*
	rm -f correctness*
	@echo "Clean complete."
