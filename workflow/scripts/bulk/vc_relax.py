import os
import sys
import adlib.bulk.vcrelax


bulk_dir = sys.argv[1]
# expecting a path that looks like this:
# path/to/dft/metal/bulk
metal = os.path.basename(os.path.dirname(bulk_dir))

adlib.bulk.vcrelax.setup_vc_relax(bulk_dir, metal=metal, lattice_constant_guess=3.6)
adlib.bulk.vcrelax.run_vc_relax(bulk_dir)
