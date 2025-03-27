#!/bin/bash

#----------------------------------------------------
# FDM MPI job script
#----------------------------------------------------
#SBATCH -J fdm_lb           # Job name
#SBATCH -o fdm_lb.out       # Name of stdout output file
#SBATCH -e fdm_lb.err       # Name of stderr output file
#SBATCH -N 5                # Total # of nodes
#SBATCH --cpus-per-task 4   # Cores per task (for OpenMP threads)
#SBATCH -n 5                # Total # of MPI tasks
#SBATCH --ntasks-per-node=1 # 1 MPI task per node
#SBATCH -t 1:00:00          # Run time (hh:mm:ss) - 1 hour
#SBATCH --nodelist=n16-80,n16-81,n16-90,n16-91,n16-92
#SBATCH --exclusive

export OMP_NUM_THREADS=4

prun ./fdm_lb 0.001953125   # Run the program with 1/512
echo "done"                 # Write this message on the output file when finished
