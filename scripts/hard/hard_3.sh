#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_3          # Job name
#SBATCH -o fdm_3.out      # Name of stdout output file(%j expands to jobId)
#SBATCH -e fdm_3.err      # Name of stderr output file(%j expands to jobId)
#SBATCH -N 1                # Total # of nodes
#SBATCH --cpus-per-task 1   # Cores per task requested
#SBATCH -n 3                # Total # of mpi tasks
#SBATCH -t 1:00:00          # Run time (hh:mm:ss) - 1 hour
#SBATCH -p lgpu             # Partition to submit to

prun ./fdm_mpi 0.00390625   # Run the program with 1/256
echo "done"                 # Write this message on the output file when finished