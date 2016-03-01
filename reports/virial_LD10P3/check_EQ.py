


from numpy import *
import matplotlib.pyplot as plt

mean_vals = []

NAS_0400 = loadtxt('../NP0400_LD10P3_C100_SF15_RT100_longer/NP0400_LD10P3_C100.ener')[:,4]
mean_NAS_0400 = mean(NAS_0400[size(NAS_0400)/2:])/(400.*25.)
mean_vals.append(mean_NAS_0400)
NAS_0500 = loadtxt('../NP0500_LD10P3_C100_NC25_SF15_RT100_longer/NP0500_LD10P3_C100.ener')[:,4]
mean_NAS_0500 = mean(NAS_0500[size(NAS_0500)/2:])/(500.*25.)
mean_vals.append(mean_NAS_0500)

NAS_0600 = loadtxt('../NP0600_LD10P3_C100_NC25_SF15_RT100_longer/NP0600_LD10P3_C100.ener')[:,4]
mean_NAS_0600 = mean(NAS_0600[size(NAS_0600)/2:])/(600.*25.)
mean_vals.append(mean_NAS_0600)

NAS_0700 = loadtxt('../NP0700_LD10P3_C100_NC25_SF15_RT100_longer/NP0700_LD10P3_C100.ener')[:,4]
mean_NAS_0700 = mean(NAS_0700[size(NAS_0700)/2:])/(700.*25.)
mean_vals.append(mean_NAS_0700)


Nt = shape(dat_0400)[0]
t = arange(Nt)*0.001

plt.clf()
plt.ion()
plt.plot(t, dat_0400[:,4]/(400.*25.), 'b-', label = 'NP0400')
plt.plot(t, dat_0500[:,4]/(500.*25.), 'r-', label = 'NP0500')
plt.plot(t, dat_0600[:,4]/(600.*25.), 'g-', label = 'NP0600')
plt.plot(t, dat_0700[:,4]/(700.*25.), 'k-', label = 'NP0700')
plt.grid()
plt.legend(loc = 'lower right')
plt.xlabel('topological dimensionless time')
plt.ylabel('fraction of NAS (%)')
mean_vals = asarray(mean_vals)
plt.yticks(mean_vals, ['%.1f' % frac_NAS for frac_NAS in mean_vals*100.])

plt.show()
         
