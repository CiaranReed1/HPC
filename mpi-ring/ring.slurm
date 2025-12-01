#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --job-name="mpi-ring"
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -t 01:00:00
#SBATCH -p shared

module load gcc
module load openmpi

echo "iterations,size,mean_time,score" > ring.dat
mpicc -O2 ring.c -o ring
for (( n = 1 ; n <= 64 ; n+=1 )); do
mpirun -n $n ./ring 10000 >> ring.dat
done 
rm ring
