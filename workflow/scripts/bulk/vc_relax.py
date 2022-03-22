import os
import sys
import adlib.bulk.vcrelax


bulk_dir = sys.argv[1]
# expecting a path that looks like this:
# path/to/dft/metal/bulk
metal = os.path.basename(os.path.dirname(bulk_dir))

nproc=16
# TODO get a reasonable guess for max_cpus from environment
adlib.bulk.vcrelax.setup_vc_relax(bulk_dir, metal=metal, lattice_constant_guess=3.6, nproc=nproc)
adlib.bulk.vcrelax.run_vc_relax(bulk_dir)
