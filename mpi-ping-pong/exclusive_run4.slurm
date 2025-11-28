#!/bin/bash

# 1 node
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --job-name="mpi-ping-pong-exclusive-4nodes"
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -t 24:00:00
#SBATCH -p multi

module load gcc
module load openmpi


mpicc -O2 variability-exclusive.c -o variability-exclusive -lm
mpirun -n 512 ./variability-exclusive > exclusive-fournode.dat
