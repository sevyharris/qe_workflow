#!/bin/bash
#SBATCH --job-name=Cu_bulk
#SBATCH --error=pw.error
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --constraint=cascadelake
#SBATCH --mincpus=32
#SbATCH --mem=100Gb


module load gcc/10.1.0
module load openmpi/4.0.5-skylake-gcc10.1
module load scalapack/2.1.0-skylake

python relax_slab.py

