#!/bin/bash

# 1 node
#SBATCH -N 1              # 1 node
#SBATCH --ntasks=1        # 1 MPI task
#SBATCH --cpus-per-task=32 # Allocate 32 cores for this job
#SBATCH --job-name="openmp-heat-color"
#SBATCH -o heat-color.out
#SBATCH -e heat-color.err
#SBATCH -t 01:00:00
#SBATCH -p shared

module load gcc
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE
echo "threads,N,time,iter" > heat-color-1024.dat
echo "threads,N,time,iter" > heat-color-1048576.dat
gcc -O2 -lm -fopenmp heat1D-pll-color.c -o r
for (( n = 1 ; n <= 32; n+=1 )); do
	./r ${n} 1024 >> heat-color-1024.dat
	./r ${n} 1048576 >> heat-color-1048576.dat
done 
rm r

