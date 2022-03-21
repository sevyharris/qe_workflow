import os
import sys
import shutil
import numpy as np

import adlib.system.calc


system_dir = sys.argv[1]
os.makedirs(system_dir, exist_ok=True)
if len(sys.argv) < 4:
    print('Assuming the slab.pwo and adsorbate.pwo files were already copied to the adsorbate-system directory')
else:
    slab_pwo = sys.argv[2]
    adsorbate_pwo = sys.argv[3]
    shutil.copy(slab_pwo, os.path.join(system_dir, 'slab.pwo'))
    shutil.copy(adsorbate_pwo, os.path.join(system_dir, 'adsorbate.pwo'))

adlib.system.calc.make_relax_script(system_dir, nproc=48)
adlib.system.calc.make_run_relax_script(system_dir, nproc=48)
adlib.system.calc.run_relax_system(system_dir)
