#!/bin/bash -l
#
# Batch script for bash users
#
#SBATCH -N 2 #the number of nodes we are using
#SBATCH -n 16 #the number of mpi ranks (1 per core), split across the nodes allocated
#SBATCH --mem=64G #up to 250G on shared queue
#SBATCH --job-name="<jobname>"
#SBATCH -o <outputfile.out> #can use /dev/null to suppress output
#SBATCH -e <outputfile.err> #can use /dev/null to suppress error output
#SBATCH -t 01:00:00 #allowed time
#SBATCH -p shared  #Queue: shared, multi, long, bigmem, test

# specify the modules you compiled the code with below
module purge #unload all modules
module load slurm/current
module load gcc
module load openmpi
module list #write a list of used modules to the outputfile

# Run the program
#add the correct command to run the program below

echo "<outputheaders>" > "<outputfile.dat>"
mpicc -O2 -lm -o "<source.c>" r
mpirun --bind-to none -n "$i" ./r >> "<outputfile.dat>"
rm r

#print some info at the end
echo "Job done, info follows..."
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode
exit
