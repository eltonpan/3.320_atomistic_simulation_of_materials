import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.pwscf import *

os.environ["PWSCF_COMMAND"] = '/home/gridsan/epan1/3320_atomistic_shared/qe-7.2/bin/pw.x'
run_qe_pwscf_from_input_file(file_path = '1C_NM.in', ncpu=8)