import os
import sys
import shutil
import numpy as np

import adlib.system.calc


# Expected arguments:
# 1. input_geometry_file - pwo file to rerun scf calculatrion for
# 2. ecutwfc in eV

# Optional arguments, but order assumed:
# 3. job_name
# 4. 

input_geometry_file = sys.argv[1]
ecutwfc = sys.argv[2]

calc_dir = os.path.join(os.path.dirname(input_geometry_file), f'scf_{ecutwfc}')
os.makedirs(calc_dir, exist_ok=True)
shutil.copy(input_geometry_file, os.path.join(calc_dir, 'input_geometry.pwo'))

# Pass in the job name
if len(sys.argv) > 3:
    job_name = sys.argv[3]
else:
    job_name = 'run_scf'


# Use low mixing_beta settings for certain metals (Cu, Ni)
low_mixing_beta=False
if len(sys.argv) > 4:
    if sys.argv[4] == 'low' or sys.argv[4] == 'low_mixing_beta':
        low_mixing_beta=True


adlib.system.calc.make_scf_script(calc_dir, nproc=48, ecutwfc=ecutwfc, low_mixing_beta=low_mixing_beta)
adlib.system.calc.make_run_scf_script(calc_dir, nproc=48, job_name=job_name)
adlib.system.calc.run_scf(calc_dir)
