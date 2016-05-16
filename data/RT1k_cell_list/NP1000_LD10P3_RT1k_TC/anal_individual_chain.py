
from numpy import *
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *
from scipy.linalg import norm
input_path      = 'tracking_individual_chain'
Np              = 1000
Nc              = 10
N_dimension     = 3
N_tot           = Np*Nc
init_cut        = 2000

chain_index     = int(sys.argv[1])
Rvec_history    = loadtxt('%s/Rvec_%06d.dat'%(input_path, chain_index))[init_cut:]

Nt = shape(Rvec_history)[0]
Rsca = zeros(Nt)
for i in range(Nt):
    Rsca[i] = norm(Rvec_history[i])
acf_Rsca = acf_gro(Rsca)
dat = zeros([size(acf_Rsca), 2])
dat[:,0] = arange(size(acf_Rsca))/100.
dat[:,1] = acf_Rsca
savetxt('%s/acf_Rsca_%06d.dat'%(input_path, chain_index), dat)
