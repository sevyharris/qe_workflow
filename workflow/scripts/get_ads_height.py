# should probable specify fcc111 for slab type

import sys
import os
import re
import glob
import numpy as np
from ase.io.espresso import read_espresso_out


if len(sys.argv) < 2:
    raise IndexError('Must specify run directory')
base_dir = sys.argv[1]

pwo_files = glob.glob(os.path.join(base_dir, 'run*', 'espresso.pwo'))

completed_files = []
for pwo_file in pwo_files:
    with open(pwo_file, 'r') as f:
        for line in reversed(f.readlines()):
            if 'JOB DONE' in line:
                completed_files.append(pwo_file)
                break


energies = np.zeros(len(completed_files))
heights = np.zeros(len(completed_files))
for i, completed_file in enumerate(completed_files):
    # get the height
    calc_file = os.path.join(os.path.dirname(completed_file), 'calc.py')
    height = 0
    with open(calc_file, 'r') as f:
        for line in f.readlines():
            m1 = re.search('height = (.*)', line)
            if m1 is not None:
                height = float(m1[1])
                break
    heights[i] = height
    with open(completed_file, 'r') as f:
        traj = list(read_espresso_out(f, index=slice(None)))
        system = traj[-1]
        energy = system.get_potential_energy()
    energies[i] = energy

best_index = np.argmin(energies)
best_height = heights[best_index]
if len(completed_files) < len(pwo_files):
    print(f'{len(completed_files)}/{len(pwo_files)} jobs completed')
    print(heights)
    print(energies)
    print(f'The best height is {best_height} with an energy of {energies[best_index]}')
else:
    print(best_height)

