#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name="mpi-gather-all"
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -t 01:00:00
#SBATCH -p shared

module load gcc
module load openmpi


mpicc -O2 gather-from-all.c -o gather-from-all 
mpirun -n 64 ./variability-exclusive 2000 256000 > gather-all.dat
mpirun -n 64 ./variability-exclusive 2000 512000 > gather-all.dat
mpirun -n 64 ./variability-exclusive 2000 1024000 > gather-all.dat
