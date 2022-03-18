import os
import sys
import numpy as np

from ase.io.espresso import read_espresso_out

import adlib.bulk.eos
import adlib.slab.calc


slab_dir = sys.argv[1]
# expecting a path that looks like this:
# path/to/dft/metal/slab
# path/to/dft/metal/bulk


if len(sys.argv) < 3:
    # load the lattice constant from eos_fine
    bulk_dir = os.path.join(os.path.dirname(slab_dir), 'bulk')
    eos_fine_dir = os.path.join(bulk_dir, 'eos_fine')
    lattice_constant = adlib.bulk.eos.analyze_eos(eos_fine_dir)
else:
    lattice_constant = float(sys.argv[2])


print(f'Using lattice constant: {lattice_constant}')
metal = os.path.basename(os.path.dirname(slab_dir))

adlib.slab.calc.make_relax_script(slab_dir, lattice_constant, metal=metal)
adlib.slab.calc.make_run_relax_script(slab_dir)
adlib.slab.calc.run_relax_slab(slab_dir)
