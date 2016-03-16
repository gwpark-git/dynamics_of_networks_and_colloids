
# from numpy import *
# import matplotlib.pyplot as plt
# from acf_fcn import *

# index_st = 2000
# dat_0400 = loadtxt('RF_NP0400_NC25_longer_2.dat')[index_st:,1]

# K_arr = [1, 2, 5, 10, 20, 50, 100, 200]
# acf_col = []
# for k in K_arr:
#     acf_col.append(acf_gro_BAV(dat_0400, k))


# colP = ['blue', 'cyan', 'green', 'brown', 'purple', 'yellow', 'red', 'black']
# lP = ['-', '-', '-', '-', '-', '-', '-', '-']
# set_alpha = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# plt.clf()
# plt.figure(figsize=(11,6))
# plt.ion()
# # for i in range(size(Nb_arr)):
# for i in range(shape(acf_col)[0]):
#     plt.plot(acf_col[i]/acf_col[i][0], lP[i], color=colP[i], label='K=%d, N_dat = %d'%(K_arr[i], 98000/K_arr[i]), alpha = set_alpha[i])

# acf_0400 = acf_col[0]
acf_0400_short = loadtxt('RF_NP0400_NC25_longer.dat')[index_st:,1]
from scipy.stats import linregress
a, y, r, p, std = linregress(range(50), log(acf_0400[:50]/acf_0400[0]))

plt.clf()
plt.ion()
plt.plot(acf_0400/acf_0400[0], 'b-', label = 'tc = %3.2f'%(-1./a))

ref_zero = asarray([[0, 0], [100, 0]])
plt.semilogy(ref_zero[:,0], ref_zero[:,1], 'k--')

plt.grid()
plt.axis([0, 200, -0.2, 1.0])
plt.legend(loc = 'upper right', numpoints=1)
plt.xlabel('topological dimensionless time')
plt.ylabel('autocorrelation function (block average method)')
plt.show()


# from scipy.stats import linregress
# a_400, y_400, r_400, p_400, std_400 = linregress(range(10), log(acf_0400[:10]/acf_0400[0]))
# a_500, y_500, r_500, p_500, std_500 = linregress(range(10), log(acf_0500[:10]/acf_0500[0]))
# a_600, y_600, r_600, p_600, std_600 = linregress(range(10), log(acf_0600[:10]/acf_0600[0]))
# a_700, y_700, r_700, p_700, std_700 = linregress(range(10), log(acf_0700[:10]/acf_0700[0]))


# plt.clf()
# plt.ion()
# plt.figure(figsize=(11,6))
# plt.plot(acf_0400/acf_0400[0], 'b-', label = 'NP0400, tc=%4.3f (based on time index)'%(-1/a_400))
# plt.plot(acf_0500/acf_0500[0], 'r-', label = 'NP0500, tc=%4.3f (based on time index)'%(-1/a_500))
# plt.plot(acf_0600/acf_0600[0], 'g-', label = 'NP0600, tc=%4.3f (based on time index)'%(-1/a_600))
# plt.plot(acf_0700/acf_0700[0], 'k-', label = 'NP0700, tc=%4.3f (based on time index)'%(-1/a_700))

# ref_zero = asarray([[0, 0], [100, 0]])
# plt.plot(ref_zero[:,0], ref_zero[:,1], 'k--')


# # plt.plot(t, acf_0400_wot/acf_0400_wot[0], 'b.-', label = 'NP0400_wot')
# # plt.plot(t, acf_0500_wot/acf_0500_wot[0], 'r-', label = 'NP0500_wot')
# # plt.plot(t, acf_0600_wot/acf_0600_wot[0], 'g-', label = 'NP0600_wot')
# # plt.plot(t, acf_0700_wot/acf_0700_wot[0], 'k-', label = 'NP0700_wot')

# plt.legend(loc = 'upper right')
# plt.grid()
# plt.ylabel('autocorrelation function')
# plt.xlabel('time index (1 = topological update)')
# plt.show()
