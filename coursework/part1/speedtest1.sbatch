#!/bin/bash

# Request resources:
#SBATCH -N 1		# number of compute nodes. 
#SBATCH --mem=10G	# memory required, up to 250G on standard nodes
#SBATCH --output=job_%j.out
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
./part1 >> speedtest_1.dat

if diff part1_${OMP_NUM_THREADS}_cores.dat serial.dat > correctness_check.txt; then
    echo "Output files match!"
    rm part1_${OMP_NUM_THREADS}_cores.dat
    rm correctness_check.txt
else
    echo "Output files differ! See correctness_check.txt for details."
fi
rm part1 
