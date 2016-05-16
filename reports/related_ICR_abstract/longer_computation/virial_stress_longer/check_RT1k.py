
from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *


# N_st_cut = 1000
# RF = loadtxt('NP0400_first_run.dat')[N_st_cut:,1]
# t_RF = arange(size(RF))/100.
# acf_RF = acf_gro(RF)
# t_acf_RF = arange(size(acf_RF))/100.
# dat = loadtxt('../NP0400_LD10P3_RT1k/NP0400_LD10P3_C100.ener')[:,4]
# t_dat = arange(size(dat))/100.

RF_RT1 = loadtxt('../../equilibrium_analysis_dt0_01/RF_NP0400_NC25_RT100.dat')

plt.close()
plt.ion()
# plt.plot(t_RF, RF, 'b-')
plt.loglog(t_acf_RF, acf_RF/acf_RF[0], 'b-')
plt.grid()
plt.show()

