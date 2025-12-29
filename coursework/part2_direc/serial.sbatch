#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -p shared
#SBATCH --time=24:00:00

SRC_BASENAME=$(basename "$BIN")
SRC_BASENAME="${SRC_BASENAME%.*}"

OUTFILE="outputs/${SRC_BASENAME}_job-${SLURM_JOB_ID}_problemsize-${PROBLEM_SIZE}.out"
mkdir -p outputs
exec > "$OUTFILE" 2>&1

echo "Running serial binary $SRC_BASENAME with problem size $PROBLEM_SIZE"

# run the program
"$BIN" "$PROBLEM_SIZE" 
