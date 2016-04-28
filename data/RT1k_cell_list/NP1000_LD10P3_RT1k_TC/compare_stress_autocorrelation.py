

from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *

# txy_virial = loadtxt('RF_NP1000_LD10P3_C100.dat')[2000:,1]
# acf_txy_virial = acf_gro(txy_virial)
# t_virial = arange(size(acf_txy_virial))/100.

# Np = 1000
# Nc = 10
# N_tot = Nc*Np
# t_RxRy = arange(50*100)/100.
# acf_av_RxRy = zeros(size(t_RxRy))
# for k in range(N_tot):
#     acf_RxRy_k = loadtxt('tracking_individual_chain/acf_RxRy_%06d.dat'%k)[:,1]
#     acf_av_RxRy[:size(acf_RxRy_k)] += acf_RxRy_k
# N_current = 8943
# acf_av_RxRy /= float(N_tot)


plt.close()
plt.ion()
plt.figure(figsize=(6,6))
plt.loglog(t_virial, acf_txy_virial/acf_txy_virial[0], 'b-', label = 'from virial')
plt.loglog(t_RxRy, acf_av_RxRy/acf_av_RxRy[0], 'r-', label = 'from individual chain')
plt.axis([-1, 2, -0.1, 1.1])
plt.legend(loc = 'upper right')
plt.xlabel(r'dimensionless time, $\beta_0 t$')
plt.ylabel('stress autocorrelation')
plt.grid()
plt.show()
