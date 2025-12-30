#!/bin/bash 
#SBATCH --partition=B1                # select partition (A1, A2, A3, A4, B1, B2, or B3) 
#SBATCH --time=24:00:00               # set time limit in HH:MM:SS 
#SBATCH --nodes=1                     # number of nodes 
#SBATCH --ntasks-per-node=20          # number of processes per node (for MPI) 
#SBATCH --cpus-per-task=1             # OMP_NUM_THREADS (for openMP) 
#SBATCH --job-name=TEST                 # job name #SBATCH --output="error.%x"
#SBATCH --output=test.out
#SBATCH --qos=short					# option: short(3days) medium(7days) long(14days)



/home/dinnercha/local/MuNRG/mybin/run_Test.sh $MCR_ROOT
