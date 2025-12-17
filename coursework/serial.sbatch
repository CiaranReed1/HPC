#!/bin/bash

# Request resources:
#SBATCH -N 1		# number of compute nodes.
#SBATCH -c 1		# number of CPU cores, one per thread, up to 128
#SBATCH --mem=1G	# memory required, up to 250G on standard nodes
#SBATCH --time=0:15:0	# time limit for job (format:  days-hours:minutes:seconds)

# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

set -e
echo "This should be run on Hamilton!"

# Modules necessary for job:
module purge
module load gcc/14.2

# clear existing text outputs
if [ -f stats.dat ]; then
	 rm stats.dat
fi
if [ -f serial ]; then
	 rm serial
fi

# compile serial.c into serial
cc -O2 serial.c -Wall -march=native -o serial -lm

# run serial
./serial
