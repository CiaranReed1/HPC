#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=32
#SBATCH --job-name="reduction-weak-scaling"
#SBATCH -o reduction-weak-scaling.%J.out
#SBATCH -e reduction-weak-scaling.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc

gcc -O1 -fopenmp reduction.c -o r
for (( n = 1 ; n <= 32 ; n+=1 )); do
OMP_NUM_THREADS=${n} ./r $((1000000 * n))
done > weak_timing.dat