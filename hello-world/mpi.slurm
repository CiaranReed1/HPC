#!/bin/bash

# 1 node
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="hello-mpi"
#SBATCH -o hello-mpi.%J.out
#SBATCH -e hello-mpi.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc
module load openmpi
mpicc -o hello-mpi mpi.c
mpirun ./hello-mpi
