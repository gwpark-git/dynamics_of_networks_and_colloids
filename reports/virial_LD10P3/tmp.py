
from numpy import *
import matplotlib.pyplot as plt
from acf_fcn import *

# index_st = 2000
# dat_0400 = loadtxt('RF_NP0400_NC25_longer.dat')[index_st:,1]
# dat_0500 = loadtxt('RF_NP0500_NC25.dat')[index_st:,1]
# dat_0600 = loadtxt('RF_NP0600_NC25.dat')[index_st:,1]
# dat_0700 = loadtxt('RF_NP0700_NC25.dat')[index_st:,1]

# cut_index = -1
# acf_0400 = acf_gro(dat_0400)[:cut_index]
# acf_0500 = acf_gro(dat_0500)[:cut_index]
# acf_0600 = acf_gro(dat_0600)[:cut_index]
# acf_0700 = acf_gro(dat_0700)[:cut_index]


t = arange(cut_index)*0.001

from scipy.stats import linregress
a_400, y_400, r_400, p_400, std_400 = linregress(range(10), log(acf_0400[:10]/acf_0400[0]))
a_500, y_500, r_500, p_500, std_500 = linregress(range(10), log(acf_0500[:10]/acf_0500[0]))
a_600, y_600, r_600, p_600, std_600 = linregress(range(10), log(acf_0600[:10]/acf_0600[0]))
a_700, y_700, r_700, p_700, std_700 = linregress(range(10), log(acf_0700[:10]/acf_0700[0]))


plt.clf()
plt.ion()
plt.figure(figsize=(11,6))
plt.plot(arange(size(acf_0400))/1000., acf_0400/acf_0400[0], 'b-', label = r'NP0400, tc=%6.4f $\tau_0$'%(-1/(1000.*a_400)))
# plt.plot(acf_0500/acf_0500[0], 'r-', label = 'NP0500, tc=%4.3f (based on time index)'%(-1/a_500))
# plt.plot(arange(size(acf_0600))/1000., acf_0600/acf_0600[0], 'g-', label = r'NP0600, tc=%4.3f $\tau_0$'%(-1/(1000.*a_600)))
plt.plot(arange(size(acf_0700))/1000., acf_0700/acf_0700[0], 'r-', label = r'NP0700, tc=%6.4f $\tau_0$'%(-1/(1000.*a_700)))

ref_zero = asarray([[0, 0], [5, 0]])
plt.plot(ref_zero[:,0], ref_zero[:,1], 'k--')


# plt.plot(t, acf_0400_wot/acf_0400_wot[0], 'b.-', label = 'NP0400_wot')
# plt.plot(t, acf_0500_wot/acf_0500_wot[0], 'r-', label = 'NP0500_wot')
# plt.plot(t, acf_0600_wot/acf_0600_wot[0], 'g-', label = 'NP0600_wot')
# plt.plot(t, acf_0700_wot/acf_0700_wot[0], 'k-', label = 'NP0700_wot')

plt.legend(loc = 'upper right')
plt.grid()
plt.ylabel('autocorrelation function')
plt.xlabel('topological dimensionless time')
plt.axis([0, 5, -0.1, 1.1])
plt.show()
