import os
import sys
import socket
from time import time
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso
#from ase.calculators.socketio import SocketIOCalculator
from ase.build import bulk
import shutil


start = time()
logfile = 'ase.log'


espresso_settings = {
    'control': {
        'verbosity': 'high',
        # 'disk_io': 'none',
        'calculation': 'vc-relax',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'degauss': 0.1,
        'ecutwfc': 500,
        # 'ecutrho': 500,
    },
    'ions': {
        'ion_dynamics': 'bfgs',
    },
    'cell': {
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'press_conv_thr': 0.005,
    }
}


cu_bulk = bulk('Cu', crystalstructure='fcc', a=3.6, cubic=True)


#unixsocket = 'ase_espresso'
#hostname = socket.gethostname()
#port = 31415  # the default port
#port = 30141  # pick a port
pw_executable = os.environ['PW_EXECUTABLE']

pseudopotentials = {
    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
    'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
    'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
    'H': 'H.pbe-kjpaw_psl.1.0.0.UPF',
}

#command = f'mpirun -np 16 {pw_executable} -in PREFIX.pwi --ipi {hostname}:{port} -nk 4 > PREFIX.pwo'
#command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi --ipi localhost:31415 -nk 4 > PREFIX.pwo'
#command = f'{pw_executable} -in PREFIX.pwi --ipi localhost:31415 > PREFIX.pwo'
#command = f'{pw_executable} -in PREFIX.pwi --ipi {unixsocket}:UNIX > PREFIX.pwo'
# command=f'aprun -n 10 -N 2 -cc depth -d 32 -j 1 {pw_executable}' -in PREFIX.pwi --ipi {hostname}:{port} -nk 5 > PREFIX.pwo',
#command=f'aprun -n 4 {pw_executable} -in PREFIX.pwi --ipi {hostname}:{port} -nk 4 > PREFIX.pwo',
#command=f'mpirun -n 8 {pw_executable} -in PREFIX.pwi --ipi {hostname}:{port} -nk 4 > PREFIX.pwo',
command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi -nk 16 > PREFIX.pwo'
print(command)


#aprun -n 8 -N 1 pw.x -nk 4 --ipi {host}:{port} --in PREFIX.pwi > PREFIX.out
espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    # kpts=(2, 2, 2),
    kpts=(16, 16, 16),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

#with SocketIOCalculator(espresso, log=sys.stdout) as calc:
#with SocketIOCalculator(espresso, log=sys.stdout, unixsocket=unixsocket) as calc:
#    cu_bulk.calc = calc
#    cu_bulk.get_potential_energy()

#with SocketIOCalculator(espresso, log=sys.stdout) as calc:
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

