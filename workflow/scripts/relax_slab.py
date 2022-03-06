import os
import sys
import socket
from time import time
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso
from ase.build import bulk, fcc111
from ase.io.espresso import read_espresso_out
from ase.io.trajectory import Trajectory
from ase.io.ulm import InvalidULMFileError



start = time()
logfile = 'ase.log'


# read in the results from the previous bulk cell relaxation
bulk_file = 'bulk.pwo'
with open(bulk_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    lattice_constant = traj[-1].cell[0][0]

print(f'Lattice constant = {lattice_constant}')

fmax = 0.01 
vacuum = 10.0
cu_slab = fcc111('Cu', size=(3, 3, 3), vacuum=vacuum, a=lattice_constant)

# Fix the bottom layer
bottom_layer = []
for i, pos in enumerate(cu_slab.get_positions()):
    if pos[-1] == vacuum:
        bottom_layer.append(cu_slab[i].index)

fix_bottom_layer = FixAtoms(indices=bottom_layer)
cu_slab.set_constraint(fix_bottom_layer)


# Restart if available
traj_file = 'slab.traj'
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            cu_slab = traj[-1]
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)


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



pw_executable = os.environ['PW_EXECUTABLE']
pseudopotentials = {
    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
    'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
    'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
    'H': 'H.pbe-kjpaw_psl.1.0.0.UPF',
}

command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi -nk 8 > PREFIX.pwo'
print(command)
espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(8, 8, 1),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)


cu_slab.calc = espresso
opt = BFGS(cu_slab, logfile=logfile, trajectory='slab.traj')
opt.run(fmax=fmax)



end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')


