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
SRC_BASENAME=$(basename "$SRC" .c)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OUTFILE="outputs/job_${SLURM_JOB_ID}_${OMP_NUM_THREADS}cores_$(SRC_BASENAME).out"
exec > "$OUTFILE" 2>&1 


gcc "$SRC" -lm -fopenmp -o bin/part1_${OMP_NUM_THREADS}_${SRC_BASENAME} 
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Running with $OMP_NUM_THREADS threads"
bin/part1_${OMP_NUM_THREADS}_${SRC_BASENAME} >> "$SPEEDFILE"

CORRECTNESS_FILE="correctness_check_${OMP_NUM_THREADS}_${SRC_BASENAME}.txt"
PROGRAM_OUTPUT="${OMP_NUM_THREADS}_cores_${SRC_BASENAME}.dat"


if diff "$PROGRAM_OUTPUT" serial.dat > "$CORRECTNESS_FILE" 2>&1; then
    echo "Output files match!"
    rm "$PROGRAM_OUTPUT" "$CORRECTNESS_FILE"
else
    echo "Output files differ! See $CORRECTNESS_FILE for details."
fi
rm part1_${OMP_NUM_THREADS}_${SRC_BASENAME}
