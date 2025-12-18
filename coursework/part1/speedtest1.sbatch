#!/bin/bash
# Request resources:
#SBATCH -N 1		# number of compute nodes. 
#SBATCH --mem=10G	# memory required, up to 250G on standard nodes
#SBATCH -e /dev/null
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --time=0:20:0	# time limit for job (format:  days-hours:minutes:seconds)
# Specify the tasks to run:
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Modules necessary for job:
module purge
module load gcc

# compile part.1c into part1
gcc -fopenmp -lm part1.c -o part1 -lm
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OUTFILE="job_${SLURM_JOB_ID}_${OMP_NUM_THREADS}cores.out"
exec > "$OUTFILE" 2>&1   # redirect all stdout/stderr

echo "Running with $OMP_NUM_THREADS threads"
./part1 >> speedtest_1.dat

CORRECTNESS_FILE="correctness_check_${OMP_NUM_THREADS}.txt"
PROGRAM_OUTPUT="part1_${OMP_NUM_THREADS}_cores.dat"

if diff "$PROGRAM_OUTPUT" serial.dat > "$CORRECTNESS_FILE" 2>&1; then
    echo "Output files match!"
    rm "$PROGRAM_OUTPUT" "$CORRECTNESS_FILE"
else
    echo "Output files differ! See $CORRECTNESS_FILE for details."
fi
rm part1 
