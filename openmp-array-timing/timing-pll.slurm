#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --job-name="timing-pll"
#SBATCH -o timing-pll.%J.out
#SBATCH -e timing-pll.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc

# unnecessary, since we've specified --cores=2 above
export OMP_NUM_THREADS=8
gcc -fopenmp -o timing-pll timing-pll.c
ulimit -s unlimited 
./timing-pll 64 50
