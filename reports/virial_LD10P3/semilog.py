
from numpy import *
import matplotlib.pyplot as plt
from acf_fcn import *

index_st = 2000
dat_0400 = loadtxt('RF_NP0400_NC25.dat')[index_st:,1]
dat_0500 = loadtxt('RF_NP0500_NC25.dat')[index_st:,1]
dat_0600 = loadtxt('RF_NP0600_NC25.dat')[index_st:,1]
dat_0700 = loadtxt('RF_NP0700_NC25.dat')[index_st:,1]

cut_index = 100
acf_0400 = acf_gro(dat_0400)[:cut_index]
acf_0500 = acf_gro(dat_0500)[:cut_index]
acf_0600 = acf_gro(dat_0600)[:cut_index]
acf_0700 = acf_gro(dat_0700)[:cut_index]

# # N_blocks = 40
# # acf_040

t = arange(cut_index)*0.001

plt.clf()
plt.ion()
plt.figure(figsize=(11,6))
plt.semilogy(acf_0400/acf_0400[0], 'b-', label = 'NP0400')
plt.semilogy(acf_0500/acf_0500[0], 'r-', label = 'NP0500')
plt.semilogy(acf_0600/acf_0600[0], 'g-', label = 'NP0600')
plt.semilogy(acf_0700/acf_0700[0], 'k-', label = 'NP0700')

ref_zero = asarray([[0, 0], [100, 0]])
plt.plot(ref_zero[:,0], ref_zero[:,1], 'k--')

# plt.plot(t, acf_0400_wot/acf_0400_wot[0], 'b.-', label = 'NP0400_wot')
# plt.plot(t, acf_0500_wot/acf_0500_wot[0], 'r-', label = 'NP0500_wot')
# plt.plot(t, acf_0600_wot/acf_0600_wot[0], 'g-', label = 'NP0600_wot')
# plt.plot(t, acf_0700_wot/acf_0700_wot[0], 'k-', label = 'NP0700_wot')

plt.legend(loc = 'lower left')
plt.grid()
plt.ylabel('autocorrelation function')
plt.xlabel('#written time step (for remapping, mutiplication with 0.001)')
# plt.show()
