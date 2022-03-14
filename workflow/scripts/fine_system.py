# should probable specify fcc111 for slab type

import sys
import os
from shutil import copyfile
import numpy as np
from ase.build import bulk, fcc111, add_adsorbate
from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms
from ase.io.espresso import read_espresso_out
from ase.io.trajectory import Trajectory
from ase.io.ulm import InvalidULMFileError
from time import time


def place_adsorbate_top(metal_slab, adsorbate, height=1.0, ads_index=0):
    # remove the adsorbate cell
    adsorbate.cell = [0, 0, 0]
    top_layer_z = np.max([pos[2] for pos in metal_slab.get_positions()])
    for i, pos in enumerate(metal_slab.get_positions()):
        if pos[2] == top_layer_z:
            # place the atom height above here
            ads_origin = pos + [0, 0, height]
            pos_difference = ads_origin - adsorbate.get_positions()[ads_index]
            adsorbate.translate(pos_difference)
            metal_slab += adsorbate
            print(f'Using {adsorbate[ads_index]} for binding atom')
            return
    print("Failed to place adsorbate")
    exit(-1)


# TODO accept optimal height guess
# if len(sys.argv) < 2:
#     raise IndexError('Must specify starting height for fine system run')
# height = float(sys.argv[1])
height = 4.0
fmax = 0.01
start = time()

logfile = 'ase.log'

slab_file = 'slab.pwo'
with open(slab_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
metal_slab = traj[-1]

adsorbate_file = 'adsorbate.pwo'
with open(adsorbate_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
adsorbate = traj[-1]

# Fix the bottom layer
bottom_layer_z = np.min([pos[2] for pos in metal_slab.get_positions()])
bottom_layer = []
for i, pos in enumerate(metal_slab.get_positions()):
    if pos[2] == bottom_layer_z:
        # print(metal_slab[i])
        bottom_layer.append(metal_slab[i].index)

fix_bottom_layer = FixAtoms(indices=bottom_layer)
metal_slab.set_constraint(fix_bottom_layer)


# place the adsorbate
element_priority = ['C', 'O', 'H']  # where does N fit?
bond_atom_index = -1
for element in element_priority:
    for atom in adsorbate:
        if atom.symbol == element:
            bond_atom_index = atom.index
            break
    if bond_atom_index > -1:
        break
print(f'bond atom is {adsorbate[bond_atom_index]}')

# 'ontop', 'bridge', 'fcc', 'hcp'
# https://wiki.fysik.dtu.dk/ase/ase/build/surface.html#ase.build.add_adsorbate
# add_adsorbate(metal_slab, adsorbate, height=height, position='ontop', mol_index=bond_atom_index)
# can't use the add_adsorbate function if you load the atoms from a pwo file... 
place_adsorbate_top(metal_slab, adsorbate, height=height, ads_index=bond_atom_index)


# Restart if available
traj_file = 'system.traj'
restart = False
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            metal_slab = traj[-1]
            restart = True
            with open(logfile, 'a') as f:
                f.write(f'Restarting from trajectory file {traj_file}\n')
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)

if not restart:
    write(f'initial_system.xyz', metal_slab)

espresso_settings = {
    'control': {
        'verbosity': 'high',
        'calculation': 'scf',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'degauss': 0.1,
        'ecutwfc': 100,
    },
}

pseudopotentials = {
    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
    'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
    'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
    'H': 'H.pbe-kjpaw_psl.1.0.0.UPF',
}

pw_executable = os.environ['PW_EXECUTABLE']
command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi > PREFIX.pwo'
calc = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(4, 4, 1),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

metal_slab.calc = calc
opt = BFGS(metal_slab, logfile=logfile, trajectory='system.traj')
opt.run(fmax=fmax)

final_energy = metal_slab.get_potential_energy()
print(f"Final Energy: {final_energy}")

# copy the file
copyfile('espresso.pwo', 'system.pwo')

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Final Energy: {final_energy}\n')
    f.write(f'Completed in {duration} seconds\n')