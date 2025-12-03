#!/bin/bash -l
#
# Batch script for bash users
#
#SBATCH -N 1 #the number of nodes we are using
#SBATCH -n 16 #the number of cores we are using
#SBATCH -J cal_pi  #a name for our job
#SBATCH -o /dev/null #the logfile
#SBATCH -e /dev/null  #the error file
#SBATCH -p shared #the queue we want to use
#SBATCH -t 00:30:00 ## the maximum run time
# Need to ensure the right MPI module is loaded -
# i.e. the same module which the program was compiled with.

# specify the modules you compiled the code with below
module purge #unload all modules
module load slurm/current
module load gcc
module load openmpi
module list #write a list of used modules to the outputfile

# Run the program
#add the correct command to run the program below
make 
echo "size,N,calculated_pi,error,time" > pll_ham.dat
for i in $(seq 1 1 $SLURM_NTASKS); do
    mpirun --bind-to none -n "$i" ./calcpi $((1000000*i)) >> pll_ham.dat
    mpirun --bind-to none -n 1 ./calcpi $((1000000*i)) >> pll_ham.dat
done
make clean

#print some info at the end
echo "Job done, info follows..."
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode
exit
