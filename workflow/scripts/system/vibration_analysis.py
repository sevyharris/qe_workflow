import os
import sys
import shutil
import numpy as np

import adlib.system.vibration


system_dir = sys.argv[1]
vib_dir = os.path.join(system_dir, 'vib')

# Pass in the job name
if len(sys.argv) > 2:
    job_name = 'vib_' + sys.argv[2]
else:
    job_name = 'vib'

# copy the system pwo file from the parent directory
os.makedirs(vib_dir, exist_ok=True)
shutil.copy(os.path.join(system_dir, 'espresso.pwo'), os.path.join(vib_dir, 'system.pwo'))

adlib.system.vibration.make_vib_analysis_script(vib_dir, nproc=48, ecutwfc=50)
adlib.system.vibration.make_run_vib_analysis_script(vib_dir, nproc=48, job_name=job_name)
adlib.system.vibration.run_vib_analysis(vib_dir)
