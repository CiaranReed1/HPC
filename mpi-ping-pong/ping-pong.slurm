#!/bin/bash

# 1 node
#SBATCH --nodes=2
#SBATCH --cores=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="mpi-ping-pong"
#SBATCH -o mpi-ping-pong.%J.out
#SBATCH -e mpi-ping-pong.%J.err
#SBATCH -t 00:10:00
#SBATCH -p shared

module load gcc
module load openmpi


mpicc -O2 hamil-ping-pong.c -o mpi-ping-pong
mpirun -n 2 ./mpi-ping-pong > two-node-pp.dat
