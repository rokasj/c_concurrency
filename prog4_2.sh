#!/bin/sh
#SBATCH -p short
#SBATCH -n2
#SBATCH -C beta
mpicc -o prog4 prog4.c
mpirun prog4 -249 250 0.0000000000001
mpirun prog4 -499 500 0.0000000000001
mpirun prog4 -249 250 0.00000000000001
mpirun prog4 -499 500 0.00000000000001
