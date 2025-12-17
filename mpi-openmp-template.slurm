#!/bin/bash -l
#
# Batch script for bash users
#

#SBATCH -c 2 #number of openmp threads per mpi rank
#SBATCH --ntasks-per-socket=1   # number of MPI ranks per CPU socket (recommend 1)
#SBATCH -N 2 #the number of nodes we are using
#SBATCH -n 4 #the number of mpi ranks (1 per socket as above, so double number of nodes)
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

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export OMP_DISPLAY_ENV=TRUE
export OMP_DISPLAY_AFFINITY=TRUE

# Run the program
#add the correct command to run the program below

echo "<outputheaders>" > "<outputfile.dat>"
mpicc -O2 -lm -fopenmp -o "<source.c>" r
mpirun --bind-to socket  --map-by socket -n "$SLURM_NTASKS" ./r >> "<outputfile.dat>"
rm r

#print some info at the end
echo "Job done, info follows..."
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode
exit
