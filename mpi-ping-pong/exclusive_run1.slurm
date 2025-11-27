#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --job-name="mpi-ping-pong-exclusive-1node"
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -t 01:00:00
#SBATCH -p multi

module load gcc
module load openmpi


mpicc -O2 variability-exclusive.c -o variability-exclusive -lm
mpirun -n 128 ./variability-exclusive > exclusive-onenode.dat
