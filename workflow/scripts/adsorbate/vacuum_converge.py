import os
import sys
import adlib.adsorbate.convergence


adsorbate_dir = sys.argv[1]
adsorbate_pwo = sys.argv[2]

adlib.adsorbate.convergence.setup_converge(adsorbate_dir, 'vacuum_converge', adsorbate_pwo=adsorbate_pwo)
adlib.adsorbate.convergence.run_converge(adsorbate_dir, 'vacuum_converge')
