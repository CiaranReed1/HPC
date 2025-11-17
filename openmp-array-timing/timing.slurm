#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --job-name="hello-openmp"
#SBATCH -o timing.%J.out
#SBATCH -e timing.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc

# unnecessary, since we've specified --cores=2 above
export OMP_NUM_THREADS=8
gcc -fopenmp -o timing timing.c
./timing 60 50
