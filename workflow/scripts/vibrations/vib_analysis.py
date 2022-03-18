import os
import sys
from time import time
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.io import read
from ase.vibrations import Vibrations


start = time()
logfile = 'ase.log'

xyz_file = '/work/westgroup/harris.se/espresso/qe_workflow/resources/initial_system.xyz'
system = read(xyz_file)


espresso_settings = {
    'control': {
        'verbosity': 'high',
        'calculation': 'scf',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.1,
        'ecutwfc': 50,
    },
}


# Fix the bottom two layers
bottom_layer = []
second_layer = []
fixed_indices = []
z_values = list(set([pos[2] for pos in system.get_positions()]))
z_values.sort()


for i, pos in enumerate(system.get_positions()):
    if pos[2] == z_values[0]:
        bottom_layer.append(system[i].index)
    if pos[2] == z_values[1]:
        second_layer.append(system[i].index)
fixed_indicies = bottom_layer + second_layer
fix_bottom_layers = FixAtoms(indices=fixed_indicies)
system.set_constraint(fix_bottom_layers)


pw_executable = os.environ['PW_EXECUTABLE']

pseudopotentials = {
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'N': 'N_ONCV_PBE-1.2.upf',
    'H': 'H_ONCV_PBE-1.2.upf',
    'Pt': 'Pt_ONCV_PBE-1.2.upf',
    'Pd': 'Pd_ONCV_PBE-1.2.upf',
}

command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi > PREFIX.pwo'
print(command)

espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(2, 2, 1),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

system.calc = espresso

vib = Vibrations(system)
vib.run()
vib.summary()

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')
