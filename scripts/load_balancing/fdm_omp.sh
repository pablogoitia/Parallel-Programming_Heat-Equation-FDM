#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_omp          # Job name
#SBATCH -o fdm_omp.out      # Name of stdout output file
#SBATCH -e fdm_omp.err      # Name of stderr output file
#SBATCH -N 8                # Total # of nodes
#SBATCH --cpus-per-task 4   # Cores per task requested
#SBATCH -n 8                # Total # of mpi tasks
#SBATCH -t 1:00:00          # Run time (hh:mm:ss) - 1 hour
#SBATCH --nodelist=n16-80,n16-81,n16-82,n16-83,n16-90,n16-91,n16-92,n16-93
#SBATCH --exclusive

export OMP_NUM_THREADS=4

prun ./fdm_omp 0.001953125  # Run the program with 1/512
echo "done"                 # Write this message on the output file when finished