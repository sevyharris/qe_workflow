# Script to relax adsorbate. Expects the adsorbate_dir and the directory with all the xyz files
import os
import sys

import adlib.adsorbate.calc


adsorbate_dir = sys.argv[1]
xyz_dir = sys.argv[2]

adlib.adsorbate.calc.setup_relax_adsorbate(adsorbate_dir, xyz_dir=xyz_dir, nproc=48)
adlib.adsorbate.calc.run_relax_adsorbate(adsorbate_dir)
