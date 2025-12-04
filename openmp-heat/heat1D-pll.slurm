#!/bin/bash

# 1 node
#SBATCH -N 1              # 1 node
#SBATCH --ntasks=1        # 1 MPI task
#SBATCH --cpus-per-task=32 # Allocate 32 cores for this job
#SBATCH --job-name="openmp-heat-pll-2-20"
#SBATCH -o heat-pll-2-20.out
#SBATCH -e heat-pll-2-20.err
#SBATCH -t 01:00:00
#SBATCH -p shared

module load gcc
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE
echo "threads,N,time,iter" > heat-pll-2-20.dat
gcc -O2 -lm -fopenmp heat1D-pll.c -o r
for (( n = 1 ; n <= 8; n+=1 )); do
	./r ${n} 1048576 >> heat-pll-2-20.dat
done 
rm r

