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

# compile part.1c into part1
SRC_BASENAME=$(basename "$SRC" .c)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
BIN="bin/part1_${OMP_NUM_THREADS}_${SRC_BASENAME}"
OUTFILE="outputs/job_${SLURM_JOB_ID}_${OMP_NUM_THREADS}cores_${SRC_BASENAME}.out"

exec > "$OUTFILE" 2>&1 


gcc "src/$SRC" -lm -fopenmp -o bin/part1_${OMP_NUM_THREADS}_${SRC_BASENAME} 
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Running $SRC with $OMP_NUM_THREADS threads"
echo "Binary: $BIN"
echo "Appending results to: $SPEEDFILE"

for run in $(seq 1 8); do
    echo "Run $run / 8"

    PROGRAM_OUTPUT="${OMP_NUM_THREADS}_cores_${SRC_BASENAME}.dat"
    CORRECTNESS_FILE="correctness_${OMP_NUM_THREADS}_${SRC_BASENAME}_run${run}.txt"

    # run program (assumes it writes PROGRAM_OUTPUT internally
    # OR writes to stdout; adapt if needed)
        (
    flock -x 200
    "$BIN" >> "$SPEEDFILE"
    ) 200>>"$SPEEDFILE.lock"

    # ---- correctness check ----
    if diff "$PROGRAM_OUTPUT" data/serial.dat > "$CORRECTNESS_FILE" 2>&1; then
        echo "Run $run: output matches serial."
        rm "$PROGRAM_OUTPUT" "$CORRECTNESS_FILE"
    else
        echo "Run $run: output differs from serial!"
        echo "See $CORRECTNESS_FILE"
    fi
done

# ---- cleanup ----
rm "$BIN"
