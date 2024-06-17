
from stereo.tools import generate_loom
import os
import sys

gem_file=sys.argv[1]
gtf_path=sys.argv[2]
out_dir =sys.argv[3]
out_file_path_loom = generate_loom(gem_path=gem_file,gtf_path=gtf_path, bin_type='bins', bin_size=40, out_dir=out_dir)
