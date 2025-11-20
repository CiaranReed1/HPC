#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=32
#SBATCH --job-name="reduction-weak-scaling-close"
#SBATCH -o reduction-weak-scaling-close.%J.out
#SBATCH -e reduction-weak-scaling-close.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc

gcc -O1 -fopenmp reduction.c -o r
for (( n = 1 ; n <= 32 ; n+=1 )); do
OMP_PLACES=cores OMP_PROC_BIND=close OMP_NUM_THREADS=${n} ./r $((1000000 * n))
done > weak_timing_close.dat