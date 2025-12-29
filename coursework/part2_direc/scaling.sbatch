#!/bin/bash
# Request resources:
#SBATCH -N 1		# number of compute nodes. 
#SBATCH --mem=10G	# memory required, up to 250G on standard nodes
#SBATCH -e /dev/null
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --time=4:0:0	# time limit for job (format:  days-hours:minutes:seconds)
# Specify the tasks to run:
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Modules necessary for job:
module purge
module load gcc
module load openmpi
# set variables
SRC_BASENAME=$(basename "$BIN")       # => 'test' (filename only)
SRC_BASENAME="${SRC_BASENAME%.*}"     # strip extension if any

OUTFILE="outputs/${SRC_BASENAME}_job-${SLURM_JOB_ID}_${SIZE}-ranks_problemsize-${PROBLEM_SIZE}.out"

#set std out of script output
exec > "$OUTFILE" 2>&1 

echo "Running $SRC_BASENAME with problem size $PROBLEM_SIZE and size $SIZE"
echo "Binary: $BIN"
echo "Appending results to: $SPEEDFILE"

#runs program
RUN_OUTPUT=$(mpirun "$BIN" "$PROBLEM_SIZE" "$RUNS")
(
    flock -x 300
    echo "$RUN_OUTPUT" >> "$SPEEDFILE"
) 300>>"$SPEEDFILE.lock"

# ---- correctness check ----
for run in $(seq 1 "$RUNS"); do
    echo "checking correctness of run $run ..."
    PROGRAM_OUTPUT="programData/test_problemsize-${PROBLEM_SIZE}_${SIZE}-ranks_run-${run}.dat"
    CORRECTNESS_FILE="correctness_${SRC_BASENAME}_${SIZE}ranks__run${run}.txt"
    if diff "$PROGRAM_OUTPUT" programData/serial.dat > "$CORRECTNESS_FILE" 2>&1; then
    echo "Run $run: output matches serial."
    rm "$PROGRAM_OUTPUT" "$CORRECTNESS_FILE"
else
    echo "Run $run: output differs from serial!"
    echo "See $CORRECTNESS_FILE"
fi
done


