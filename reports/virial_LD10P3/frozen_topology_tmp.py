
from numpy import *
import matplotlib.pyplot as plt
dat_top = loadtxt('RF_NP0700_NC25.dat')[:,1]
t_top = arange(size(dat_top))*0.001
dat_frozen = loadtxt('RF_frozen_NP0700_NC25.dat')[:,1]
t_frozen = t_top[-1] + arange(size(dat_frozen))*0.1*0.001

plt.clf()
plt.ion()
plt.plot(t_top, dat_top, 'b-', label = 'topological update (NP=700)')
plt.plot(t_frozen, dat_frozen, 'r-', label = 'frozen topology (NP=700)')
plt.grid()
plt.xlabel('topological dimensionless time')
plt.ylabel('shear stress function')
plt.legend(loc = 'upper right')
plt.show()
