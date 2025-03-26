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
for n in 0.015625 0.012401571 0.01083377 0.009843133 0.009137555 0.008598769 0.008168093 0.0078125 0.007511717 0.007252483 0.007025692 0.006824847 0.006645162 0.00648302 0.006335627 0.006200785
do
    echo "Running with value: $n"
    srun ./fdm_sec $n
done

echo "done" # Write this message on the output file when finished