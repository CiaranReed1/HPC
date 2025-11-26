#!/bin/bash

# 1 node
#SBATCH --nodes=1
#SBATCH --cores=2
#SBATCH --job-name="mpi-ping-pong"
#SBATCH -o mpi-ping-pong.%J.out
#SBATCH -e mpi-ping-pong.%J.err
#SBATCH -t 00:01:00
#SBATCH -p shared

module load gcc
module load mpi

mpicc -O2 hamil-ping-pong.c -o mpi-ping-pong
mpirun -n 2 ./mpi-ping-pong > mpi-ping-pong.dat