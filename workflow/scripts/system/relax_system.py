import os
import sys
import shutil
import numpy as np

import adlib.system.calc


# Expected arguments:
# 1. system_dir - path to the system directory where the calculation will take place
# 2. slab_pwo
# 3. adsorbate_pwo

# Optional arguments:
# 4. job_name
# 5. low_mixing_beta

system_dir = sys.argv[1]

os.makedirs(system_dir, exist_ok=True)
if len(sys.argv) < 4:
    print('Assuming the slab.pwo and adsorbate.pwo files were already copied to the adsorbate-system directory')
else:
    slab_pwo = sys.argv[2]
    adsorbate_pwo = sys.argv[3]
    shutil.copy(slab_pwo, os.path.join(system_dir, 'slab.pwo'))
    shutil.copy(adsorbate_pwo, os.path.join(system_dir, 'adsorbate.pwo'))

# Pass in the job name
if len(sys.argv) > 4:
    job_name = sys.argv[4]
else:
    job_name = 'relax_system'

# Use low mixing_beta settings for certain metals (Cu, Ni)
low_mixing_beta=False
if len(sys.argv) > 5:
    if sys.argv[5] == 'low' or sys.argv[5] == 'low_mixing_beta':
        low_mixing_beta=True


adlib.system.calc.make_relax_script(system_dir, nproc=48, ecutwfc=50, low_mixing_beta=low_mixing_beta)
adlib.system.calc.make_run_relax_script(system_dir, nproc=48, job_name=job_name)
adlib.system.calc.run_relax_system(system_dir)
