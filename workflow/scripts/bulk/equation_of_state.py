# Script to set up energy calculations for an Energy vs. Lattice constant plot
import os
import sys
from time import time

import numpy as np
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.build import bulk

import shutil
# import job_manager
import run_bulk_energy


# BASE_DIR = sys.argv[1]
BASE_DIR = '/work/westgroup/harris.se/espresso/qe_workflow/results/bulk/equation_of_state/'
lattice_constants = np.linspace(3.5, 3.7, 21)

for i, a in enumerate(lattice_constants):
    run_bulk_energy.setup_energy_calc(os.path.join(BASE_DIR, f'run{i}'), a, array_job=True)

run_bulk_energy.make_scf_run_file_array(BASE_DIR, i)
