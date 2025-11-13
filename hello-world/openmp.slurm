#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --job-name="hello-openmp"
#SBATCH -o hello-openmp.%J.out
#SBATCH -e hello-openmp.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc

# unnecessary, since we've specified --cores=2 above
export OMP_NUM_THREADS=8
gcc -fopenmp -o hello-openmp openmp.c
./hello-openmp
