import os
import sys
from time import time
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.build import bulk
import shutil


start = time()
logfile = 'ase.log'


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
        'ecutwfc': 500,
    },
}


cu_bulk = bulk('Cu', crystalstructure='fcc', a=3.6, cubic=True)

pw_executable = os.environ['PW_EXECUTABLE']

use_oncv = True

if use_oncv:
    pseudopotentials = {
        'C': 'C_ONCV_PBE-1.2.upf',
        'Cu': 'Cu_ONCV_PBE-1.2.upf',
        'O': 'O_ONCV_PBE-1.2.upf',
        'N': 'N_ONCV_PBE-1.2.upf',
        'H': 'H_ONCV_PBE-1.2.upf',
    }
else:
    pseudopotentials = {
        'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
        'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
        'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
        'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
        'H': 'H.pbe-kjpaw_psl.1.0.0.UPF',
    }

command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi -nk 16 > PREFIX.pwo'
print(command)

espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(16, 16, 16),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

cu_bulk.calc = espresso
cu_bulk.get_potential_energy()


end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')


# if it worked, copy the file
job_done = False
with open('espresso.pwo', 'r') as f:
    for line in f.readlines():
        if 'JOB DONE' in line:
            job_done = True          
            break

if job_done:
    shutil.copyfile('espresso.pwo', 'bulk.pwo')
    print('JOB DONE, copying to bulk.pwo')
else:
    print('incomplete, not copying')

