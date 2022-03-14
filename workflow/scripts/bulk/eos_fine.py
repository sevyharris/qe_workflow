import os
import sys
import numpy as np

from ase.io.espresso import read_espresso_out

import adlib.bulk.eos


bulk_dir = sys.argv[1]

if len(sys.argv) < 3:
    eos_coarse_dir = os.path.join(bulk_dir, 'eos_coarse')
    lattice_constant_guess = adlib.bulk.eos.analyze_eos(eos_coarse_dir)
else:
    lattice_constant_guess = float(sys.argv[2])


# expecting a path that looks like this:
# path/to/dft/metal/bulk
metal = os.path.basename(os.path.dirname(bulk_dir))

adlib.bulk.eos.setup_eos_fine(
    bulk_dir,
    metal=metal,
    lattice_constant_guess=lattice_constant_guess,
)
calc_dir = os.path.join(bulk_dir, 'eos_fine')
adlib.bulk.eos.run_eos(calc_dir)
