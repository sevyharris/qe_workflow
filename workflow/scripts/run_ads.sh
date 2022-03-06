#!/bin/bash
#SBATCH --error=pw.error
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --mincpus=4
#SBATCH --constraint=cascadelake


module load gcc/10.1.0
module load openmpi/4.0.5-skylake-gcc10.1
module load scalapack/2.1.0-skylake

python relax_adsorbate.py

