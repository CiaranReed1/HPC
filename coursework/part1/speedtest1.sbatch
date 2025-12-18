#!/bin/bash

# Request resources:
#SBATCH -N 1		# number of compute nodes. 
#SBATCH -c 1		# number of CPU cores, one per thread, up to 128
#SBATCH --mem=1G	# memory required, up to 250G on standard nodes
#SBATCH --output=job_%j.out
#SBATCH --time=0:15:0	# time limit for job (format:  days-hours:minutes:seconds)

# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Modules necessary for job:
module purge
module load gcc

# compile part.1c into part1
gcc -fopenmp -lm part1.c -o part1 -lm

# run part1 with 1 thread
echo "cores,time_seconds" > speedtest_1.dat

export OMP_NUM_THREADS=1
./part1 >> speedtest_1.dat

diff part1.dat serial.dat > correctness_check.txt
if diff part1.dat serial.dat > correctness_check.txt; then
    echo "Output files match!"
    rm part1.dat
    rm correctness_check.txt
else
    echo "Output files differ! See correctness_check.txt for details."
fi
rm part1 
