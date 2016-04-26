
from numpy import *
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *

input_path      = 'tracking_individual_chain'
Np              = 1000
Nc              = 10
N_dimension     = 3

N_tot           = Np*Nc

chain_index     = int(sys.argv[1])
Rvec_history    = loadtxt('%s/Rvec_%06d.dat'%(input_path, chain_index))
