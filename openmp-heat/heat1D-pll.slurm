#!/bin/bash

# 1 node
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --job-name="openmp-heat-pll"
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -t 01:00:00
#SBATCH -p shared

module load gcc
echo "threads,N,time,iter" > heat-pll.dat
gcc -O2 -lm -fopenmp openmp-heat-pll.c -o r
for (( n = 1 ; n <= 32 ; n+=1 )); do
OMP_PLACES=cores OMP_PROC_BIND=spread ./r ${n} 4096 >> heat-pll.dat
done 
rm r

