
# from numpy import *
# import matplotlib.pyplot as plt

# import sys
# sys.path.append('../../../../post_processing')
# from acf_fcn import *
# normal_factor = 2.*1000.

# sf_NP0400 = loadtxt('RF_NP0400_NC25_RT100.dat')[1000:,1]
# t_sf_NP0400 = arange(size(sf_NP0400))/100.
# acf_NP0400 = acf_gro(sf_NP0400/normal_factor)
# t_acf_NP0400 = arange(size(acf_NP0400))/100.

# sf_NP0600 = loadtxt('RF_NP0600_NC25_RT100.dat')[1000:,1]
# t_sf_NP0600 = arange(size(sf_NP0600))/100.
# acf_NP0600 = acf_gro(sf_NP0600/normal_factor)
# t_acf_NP0600 = arange(size(acf_NP0600))/100.

# sf_NP0800 = loadtxt('RF_NP0800_NC25_RT100.dat')[1000:,1]
# t_sf_NP0800 = arange(size(sf_NP0800))/100.
# acf_NP0800 = acf_gro(sf_NP0800/normal_factor)
# t_acf_NP0800 = arange(size(acf_NP0800))/100.


# sf_NP1000 = loadtxt('RF_NP1000.dat')[1000:,1]
# t_sf_NP1000 = arange(size(sf_NP1000))/100.
# normal_factor = 2.*1000.
# acf_NP1000 = acf_gro(sf_NP1000/normal_factor)
# t_acf_NP1000 = arange(size(acf_NP1000))/100.

# sf_NP1200 = loadtxt('RF/RF_NP1200_RT100_dt0_01_longer_combined.dat')[1000:,1]
# t_sf_NP1200 = arange(size(sf_NP1200))/100.
# acf_NP1200 = acf_gro(sf_NP1200/normal_factor)
# t_acf_NP1200 = arange(size(acf_NP1200))/100.

# sf_NP1400 = loadtxt('RF/RF_NP1400_RT100_dt0_01_longer_combined.dat')[1000:,1]
# t_sf_NP1400 = arange(size(sf_NP1400))/100.
# acf_NP1400 = acf_gro(sf_NP1400/normal_factor)
# t_acf_NP1400 = arange(size(acf_NP1400))/100.

colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple']
symP = ['o', 's', 'D', '^', 'v', 'h']

peak_info = []


def find_02(dat):
    for i,x in enumerate(dat):
        if x > 0.2:
            return i
    return i

def find_val(dat, val):
    for i,x in enumerate(dat/dat[0]):
        if x < val:
            return i
    return i

val = 0.2
t_acf_NP0400_reduced, acf_NP0400_reduced = t_acf_NP0400[:100], acf_NP0400[:100]
tsf_NP0400 = t_acf_NP0400[find_val(acf_NP0400_reduced, val)]
t_acf_NP0600_reduced, acf_NP0600_reduced = t_acf_NP0600[:120], acf_NP0600[:120]
tsf_NP0600 = t_acf_NP0600[find_val(acf_NP0600_reduced, val)]

t_acf_NP0800_reduced, acf_NP0800_reduced = t_acf_NP0800[:160], acf_NP0800[:160]
tsf_NP0800 = t_acf_NP0800[find_val(acf_NP0800_reduced, val)]

t_acf_NP1000_reduced, acf_NP1000_reduced = t_acf_NP1000[:430], acf_NP1000[:430]
tsf_NP1000 = t_acf_NP1000[find_val(acf_NP1000_reduced, val)]

t_acf_NP1200_reduced, acf_NP1200_reduced = t_acf_NP1200[:700], acf_NP1200[:700]
tsf_NP1200 = t_acf_NP1200[find_val(acf_NP1200_reduced, val)]

t_acf_NP1400_reduced, acf_NP1400_reduced = t_acf_NP1400[:1400], acf_NP1400[:1400]
tsf_NP1400 = t_acf_NP1400[find_val(acf_NP1400_reduced, val)]

NP = [400, 600, 800, 1000, 1200, 1400]
sup_time = asarray([tsf_NP0400, tsf_NP0600, tsf_NP0800, tsf_NP1000, tsf_NP1200, tsf_NP1400])

# plt.close()
# plt.ion()
# plt.plot(NP, sup_time, 'bo-')
# plt.show()

plt.close()
plt.ion()
plt.figure(figsize=(11,6))
plt.loglog(t_acf_NP0400_reduced/tsf_NP0400, acf_NP0400_reduced/acf_NP0400_reduced[0], '-', linewidth=2, alpha=0.7, color=colP[0], label = r'$\nu_0=10/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP0400[0])
plt.loglog(t_acf_NP0600_reduced/tsf_NP0600, acf_NP0600_reduced/acf_NP0600_reduced[0], '-', linewidth=2, alpha=0.7, color=colP[1], label = r'$\nu_0=15/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP0600[0])
plt.loglog(t_acf_NP0800_reduced/tsf_NP0800, acf_NP0800_reduced/acf_NP0800_reduced[0], '-', linewidth=2, alpha=0.7, color=colP[2], label = r'$\nu_0=20/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP0800[0])
plt.loglog(t_acf_NP1000_reduced/tsf_NP1000, acf_NP1000_reduced/acf_NP1000_reduced[0], '-', linewidth=2, alpha=0.7, color=colP[3], label = r'$\nu_0=25/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP1000[0])
plt.loglog(t_acf_NP1200_reduced/tsf_NP1200, acf_NP1200_reduced/acf_NP1200_reduced[0], '-', linewidth=2, alpha=0.7, color=colP[4], label = r'$\nu_0=30/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP1200[0])
plt.loglog(t_acf_NP1400_reduced/tsf_NP1400, acf_NP1400_reduced/acf_NP1400_reduced[0], '-', linewidth=2, alpha=0.7, color=colP[5], label = r'$\nu_0=35/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP1400[0])
plt.xlabel(r'topological time (normalized by $\tau_0$)')
plt.ylabel('dimensionless stress autocorrelation for bridges')
plt.axis([0, 20, 10**-2, 10**0])

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('x-small')
plt.grid()

plt.legend(loc = 'lower left', numpoints=1, ncol=2, prop=fontP)
plt.show()
