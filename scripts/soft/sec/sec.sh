#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_sec          # Job name
#SBATCH -o fdm_sec.out      # Name of stdout output file
#SBATCH -e fdm_sec.err      # Name of stderr output file
#SBATCH -N 1                # Total # of nodes
#SBATCH --cpus-per-task 1   # Cores per task requested
#SBATCH -t 5:00:00          # Run time (hh:mm:ss) - 5 hours
#SBATCH -p lgpu             # Partition to submit to

# Run for different values of n
for n in 0.00390625 0.003100393 0.002708442 0.002460783 0.002284389 0.002149692 0.002042023 0.001953125 0.001877929 0.001813121 0.001756423 0.001706212 0.001661291 0.001620755 0.001583907 0.001550196
do
    echo "Running with value: $n"
    srun ./fdm_sec $n
done

echo "done" # Write this message on the output file when finished