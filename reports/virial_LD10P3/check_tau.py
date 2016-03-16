
from numpy import *
import matplotlib.pyplot as plt
from acf_fcn import *

index_st = 2000
dat_0400 = loadtxt('RF_NP0400_NC25.dat')[index_st:,1]
dat_0500 = loadtxt('RF_NP0500_NC25.dat')[index_st:,1]
dat_0600 = loadtxt('RF_NP0600_NC25.dat')[index_st:,1]
dat_0700 = loadtxt('RF_NP0700_NC25.dat')[index_st:,1]

plt.clf()
plt.figure(figsize=(16,12))
plt.ion()
t = arange(8000) + 2000
plt.plot(t, dat_0400, 'b-', label = 'Np=400')
plt.plot(t, dat_0500 + 300, 'r-', label = 'Np=500, shifted by 300')
plt.plot(t, dat_0600 + 600, 'g-', label = 'Np=600, shifted by 600')
plt.plot(t, dat_0700 + 900, 'k-', label = 'Np=900, shifted by 900')
plt.grid()
plt.legend(loc = 'upper right', numpoints=1)
plt.yticks([0, 300, 600, 900])
# plt.axis([2000, 2100, -200, 1600])
plt.xlabel('time index (1 = topological update)')
plt.ylabel('shifted shear stress')
plt.show()
