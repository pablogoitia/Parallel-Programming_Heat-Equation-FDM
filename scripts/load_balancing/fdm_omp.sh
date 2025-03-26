#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_omp          # Job name
#SBATCH -o fdm_omp.out      # Name of stdout output file
#SBATCH -e fdm_omp.err      # Name of stderr output file
#SBATCH -N 4                # Total # of nodes
#SBATCH --cpus-per-task 4   # Cores per task (for OpenMP threads)
#SBATCH -n 4                # Total # of MPI tasks
#SBATCH --ntasks-per-node=1 # 1 MPI task per node
#SBATCH -t 1:00:00          # Run time (hh:mm:ss) - 1 hour
#SBATCH --nodelist=n16-82,n16-83,n16-92,n16-93
#SBATCH --exclusive


export OMP_NUM_THREADS=4

prun ./fdm_omp 0.0009765625  # Run the program with 1/1024
echo "done"                  # Write this message on the output file when finished