

# from numpy import *
# import matplotlib.pyplot as plt
# from acf_fcn import *

# dat = loadtxt('RF_frozen_NP0700_NC25.dat')[:,1]
# acf_dat = acf_gro(dat[10:])


# dat_top = loadtxt('RF_NP0700_NC25.dat')[2000:,1]
# acf_top = acf_gro(dat_top)
t_top = arange(size(acf_top))*0.001

Nt = size(acf_dat)
t = arange(Nt)*0.1*0.001


from scipy.stats import linregress
a, y, r, p, std = linregress(t[:50], log(acf_dat[:50]))
a_top, y_top, r_top, p_top, std_top = linregress(t_top[:10], log(acf_top[:10]))

plt.clf()
plt.ion()
# plt.figure(figsize=(11,6))
plt.plot(t_top, acf_top/acf_top[0], 'b-', label = r'Update topology, initial tc = %6.4f $\tau_0$'%(-1./a_top))
plt.plot(t, acf_dat/acf_dat[0],'r-', label = r'Frozen topology, initial tc = %6.4f $\tau_0$'%(-1./a))
plt.grid()
plt.xlabel('topological dimensionless time')
plt.ylabel('normalized autocorrelation')
plt.legend(loc = 'center right')
plt.axis([0, 0.5, -0.1, 1.1])
plt.show()
