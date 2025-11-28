#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=2
#SBATCH --ntasks-per-node=2
#SBATCH --job-name="mpi-ping-pong"
#SBATCH -o mpi-ping-pong.%J.out
#SBATCH -e mpi-ping-pong.%J.err
#SBATCH -t 00:10:00
#SBATCH -p shared

module load gcc
module load openmpi


mpicc -O2 variability.c -o variability -lm
mpirun -n 2 ./variability > variability-onenode.dat