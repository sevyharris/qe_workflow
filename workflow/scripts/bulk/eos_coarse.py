import os
import sys
import numpy as np

from ase.io.espresso import read_espresso_out

import adlib.bulk.eos


bulk_dir = sys.argv[1]

if len(sys.argv) < 3:
    vc_relax_file = os.path.join(bulk_dir, 'vc_relax', 'espresso.pwo')
    with open(vc_relax_file, 'r') as f:
        traj = list(read_espresso_out(f, index=slice(None)))
        atoms = traj[-1]
        lattice_constant_guess = atoms.get_distances(0, 1)[0] * np.sqrt(2)
else:
    lattice_constant_guess = float(sys.argv[2])

print(f'Lattice constant guess: {lattice_constant_guess}')
# expecting a path that looks like this:
# path/to/dft/metal/bulk
metal = os.path.basename(os.path.dirname(bulk_dir))

adlib.bulk.eos.setup_eos_coarse(
    bulk_dir,
    metal=metal,
    lattice_constant_guess=lattice_constant_guess,
)
calc_dir = os.path.join(bulk_dir, 'eos_coarse')
adlib.bulk.eos.run_eos(calc_dir)

