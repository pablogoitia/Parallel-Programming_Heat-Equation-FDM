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
for n in 4 8 16 32 64 128 256 512 1024 2048 4096
do
    value=$(echo "1/$n" | bc -l)
    echo "Running with value: $value"
    srun ./fdm_sec $value
done

echo "done" # Write this message on the output file when finished