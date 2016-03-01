
from numpy import *
import matplotlib.pyplot as plt
from acf_fcn import *

index_st = 2000
dat_0400 = loadtxt('RF_NP0400_NC25.dat')[index_st:,1]
dat_0500 = loadtxt('RF_NP0500_NC25.dat')[index_st:,1]
dat_0600 = loadtxt('RF_NP0600_NC25.dat')[index_st:,1]
dat_0700 = loadtxt('RF_NP0700_NC25.dat')[index_st:,1]

cut_index = 100
acf_0400 = acf_np(dat_0400)[:cut_index]
acf_0500 = acf_np(dat_0500)[:cut_index]
acf_0600 = acf_np(dat_0600)[:cut_index]
acf_0700 = acf_np(dat_0700)[:cut_index]

# N_blocks = 40
# acf_040

t = arange(cut_index)*0.001

plt.clf()
plt.ion()
plt.plot(acf_0400/acf_0400[0], 'b-', label = 'NP0400')
plt.plot(acf_0500/acf_0500[0], 'r-', label = 'NP0500')
plt.plot(acf_0600/acf_0600[0], 'g-', label = 'NP0600')
plt.plot(acf_0700/acf_0700[0], 'k-', label = 'NP0700')
# plt.plot(t, acf_0400_wot/acf_0400_wot[0], 'b.-', label = 'NP0400_wot')
# plt.plot(t, acf_0500_wot/acf_0500_wot[0], 'r-', label = 'NP0500_wot')
# plt.plot(t, acf_0600_wot/acf_0600_wot[0], 'g-', label = 'NP0600_wot')
# plt.plot(t, acf_0700_wot/acf_0700_wot[0], 'k-', label = 'NP0700_wot')

plt.legend(loc = 'upper right')
plt.grid()
plt.ylabel('autocorrelation function')
plt.xlabel('#written time step (for remapping, mutiplication with 0.001)')
plt.title('check without time average')
plt.show()
