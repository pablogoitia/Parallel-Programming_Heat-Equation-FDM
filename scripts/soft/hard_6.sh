#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_6          # Job name
#SBATCH -o fdm_6.out      # Name of stdout output file(%j expands to jobId)
#SBATCH -e fdm_6.err      # Name of stderr output file(%j expands to jobId)
#SBATCH -N 2                # Total # of nodes
#SBATCH --cpus-per-task 1   # Cores per task requested
#SBATCH -n 6                # Total # of mpi tasks
#SBATCH -t 1:00:00          # Run time (hh:mm:ss) - 1 hour
#SBATCH -p lgpu             # Partition to submit to

prun ./fdm_mpi 0.002149692
echo "done"                 # Write this message on the output file when finished