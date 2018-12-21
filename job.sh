#!/bin/bash

#SBATCH -J flash4-test
#SBATCH -p debug
#SBATCH --mem-per-cpu=50MB
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
##SBATCH --exclusive
#SBATCH -t 0:30:00
#SBATCH --no-requeue

##module load openmpi/2.1.1/b1
module load impi/2017
module load hdf5/1.8.17/b4
module load hypre/2.11.1/b2
module load intel/2017

mpirun -np 16 ./flash4 > flash4.log
