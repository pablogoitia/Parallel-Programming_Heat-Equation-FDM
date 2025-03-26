#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_sec          # Job name
#SBATCH -o fdm_sec.out      # Name of stdout output file(%j expands to jobId)
#SBATCH -e fdm_sec.err      # Name of stderr output file(%j expands to jobId)
#SBATCH -N 1                # Total # of nodes
#SBATCH --cpus-per-task 1   # Cores per task requested
#SBATCH -t 5:00:00          # Run time (hh:mm:ss) - 5 hours
#SBATCH -p lgpu             # Partition to submit to

# Run for different values of n (powers of 2)
# Probably it will not finish, but we will let it compute as much as it can
for n in 0.25 0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625 0.001953125 0.0009765625
do
    echo "Running with value: $n"
    srun ./fdm_sec $n
done

echo "done" # Write this message on the output file when finished