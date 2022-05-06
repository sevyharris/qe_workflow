import os
import sys
import shutil
import numpy as np

import adlib.system.vibration


system_dir = sys.argv[1]
vib_dir = os.path.join(system_dir, 'vib')

# expects
# 1. system dir
# 2. system name for job name
# 3. mixing beta

# Pass in the job name
if len(sys.argv) > 2:
    job_name = 'vib_' + sys.argv[2]
else:
    job_name = 'vib'

# Use low mixing_beta settings for certain metals (Cu, Ni)
low_mixing_beta = False
if len(sys.argv) > 3:
    if sys.argv[3] == 'low' or sys.argv[3] == 'low_mixing_beta':
        low_mixing_beta = True


# copy the system pwo file from the parent directory
os.makedirs(vib_dir, exist_ok=True)
shutil.copy(os.path.join(system_dir, 'espresso.pwo'), os.path.join(vib_dir, 'system.pwo'))

adlib.system.vibration.make_vib_analysis_script(vib_dir, nproc=16, ecutwfc=50, low_mixing_beta=low_mixing_beta)
adlib.system.vibration.make_run_vib_analysis_script(vib_dir, nproc=16, job_name=job_name)
adlib.system.vibration.run_vib_analysis(vib_dir)
