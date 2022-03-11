import os
import sys

import numpy as np


BASE_DIR = sys.argv[1]
# BASE_DIR = '/work/westgroup/harris.se/espresso/qe_workflow/results/bulk/'

# Try one parameter at a time
# ecuts = np.logspace(1, 3.5, 11)  # 10 - 3163 
# kpts = [k for k in range(1, 22)] 
# smearing = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]

ecuts = np.logspace(1, 3.5, 11)  # 10 - 3163 
kpts = [9]
smearing = [0.1]


def make_scf_run_file(dest_dir, N_runs):
    bash_filename = os.path.join(dest_dir, 'run_qe_jobs.sh')
    run_i_dir = os.path.abspath(os.path.join(dest_dir, 'run$SLURM_ARRAY_TASK_ID'))
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('#SBATCH --time=24:00:00\n')
        f.write('#SBATCH --job-name=lattice_converge\n')
        f.write('#SBATCH --mem=40Gb\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --ntasks=16\n')
        f.write('#SBATCH --partition=short,west\n')
        f.write(f'#SBATCH --array=0-{N_runs - 1}\n\n')
        f.write('module load gcc/10.1.0\n')
        f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
        f.write('module load scalapack/2.1.0-skylake\n\n')
        f.write(f'cd {run_i_dir}\n')
        f.write(f'python calc.py\n')

def make_scf_calc_file(ecutwfc, kpt, smear, calc_dir):

    python_file_lines = [
        "import os",
        "import sys",
        "from time import time",
        "from ase.calculators.espresso import Espresso",
        "from ase.build import bulk",
        "",
        "",
        "start = time()",
        "logfile = 'ase.log'",
        "",
        "",
        "espresso_settings = {",
        "    'control': {",
        "        'verbosity': 'high',",
        "        'calculation': 'scf',",
        "    },",
        "    'system': {",
        "        'input_dft': 'BEEF-VDW',",
        "        'occupations': 'smearing',",
        "        'smearing': 'mv',",
        f"        'degauss': {smear},",
        f"        'ecutwfc': {ecutwfc},",
        "    },",
        "}",
        "",
        "",
        "cu_bulk = bulk('Cu', crystalstructure='fcc', a=3.6, cubic=True)",
        "",
        "pw_executable = os.environ['PW_EXECUTABLE']",
        "",
        "pseudopotentials = {",
        "    'C': 'C_ONCV_PBE-1.2.upf',",
        "    'Cu': 'Cu_ONCV_PBE-1.2.upf',",
        "    'O': 'O_ONCV_PBE-1.2.upf',",
        "    'N': 'N_ONCV_PBE-1.2.upf',",
        "    'H': 'H_ONCV_PBE-1.2.upf',",
        "}",
        "",
        "command = f'mpirun -np 16 {pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "print(command)",
        "",
        "espresso = Espresso(",
        "    command=command,",
        "    pseudopotentials=pseudopotentials,",
        "    tstress=True,",
        "    tprnfor=True,",
        f"    kpts=({kpt}, {kpt}, {kpt}),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "cu_bulk.calc = espresso",
        "energy = cu_bulk.get_potential_energy()",
        "",
        "",
        "end = time()",
        "duration = end - start",
        "",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Energy: {energy} eV\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
        "",
    ]

    calc_filename = os.path.join(calc_dir, f'calc.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])

i = 0
for ecutwfc in ecuts:
    for kpt in kpts:
        for smear in smearing:
            calc_dir = os.path.join(BASE_DIR, f'run{i}')
            os.makedirs(calc_dir, exist_ok=True)
            make_scf_calc_file(ecutwfc, kpt, smear, calc_dir)
            i += 1


make_scf_run_file(BASE_DIR, i)
