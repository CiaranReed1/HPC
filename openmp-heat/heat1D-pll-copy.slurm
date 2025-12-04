#!/bin/bash

# 1 node
#SBATCH -N 1              # 1 node
#SBATCH --ntasks=1        # 1 MPI task
#SBATCH --cpus-per-task=32 # Allocate 32 cores for this job
#SBATCH --job-name="openmp-heat-copy"
#SBATCH -o heat-copy.out
#SBATCH -e heat-copy.err
#SBATCH -t 01:00:00
#SBATCH -p shared

module load gcc
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE
echo "threads,N,time,copy_time,iter" > heat-copy-1024.dat
echo "threads,N,time,copy_time,iter" > heat-copy-1048576.dat
gcc -O2 -lm -fopenmp heat1D-pll-copy.c -o r
for (( n = 1 ; n <= 8; n+=1 )); do
	./r ${n} 1024 >> heat-copy-1024.dat
	./r ${n} 1048576 >> heat-copy-1048576.dat
done 
rm r

