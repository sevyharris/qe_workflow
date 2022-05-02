#!/bin/bash

#SBATCH --time=unlimited
#SBATCH --job-name=vib
#SBATCH --mem=40Gb
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=32
#SBATCH --partition=west
module load gcc/10.1.0
module load openmpi/4.0.5-skylake-gcc10.1
module load scalapack/2.1.0-skylake

python  vib_analysis.py

