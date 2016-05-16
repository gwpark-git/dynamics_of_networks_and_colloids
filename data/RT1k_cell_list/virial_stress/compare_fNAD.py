
from numpy import *
import matplotlib.pyplot as plt

Nc = 10.
fac_index_to_time = 100.
frac_to_percent = 100.
fNAS_NP1000 = frac_to_percent*append(append(loadtxt('../NP1000_LD10P3_RT1k/NP1000_LD10P3_C100.ener')[:,4]/(Nc*1000.),
                                            loadtxt('../NP1000_LD10P3_RT1k_cont/NP1000_LD10P3_C100.ener')[:,4]/(Nc*1000.)),
                                     loadtxt('../NP1000_LD10P3_RT1k_cont_2/NP1000_LD10P3_C100.ener')[:,4]/(Nc*1000.))

t_NP1000 = arange(size(fNAS_NP1000))/fac_index_to_time
mean_fNAS_NP1000 = mean(fNAS_NP1000[int(10*fac_index_to_time):])

fNAS_NP1728 = frac_to_percent*append(loadtxt('../NP1728_LD12P3_RT1k/NP1728_LD12P3_C100.ener')[:,4]/(Nc*1728.),
                                     loadtxt('../NP1728_LD12P3_RT1k_cont/NP1728_LD12P3_C100.ener')[:,4]/(Nc*1728.))
t_NP1728 = arange(size(fNAS_NP1728))/fac_index_to_time
mean_fNAS_NP1728 = mean(fNAS_NP1728[int(10*fac_index_to_time):])


plt.close()
plt.figure(figsize=(11, 6))
plt.ion()
plt.plot(t_NP1000, fNAS_NP1000, 'b-', label = r'Np=1000, $\tilde{V}= 10^3$, mean(fNAS)=%4.2f (%%)'%mean_fNAS_NP1000)
plt.plot(t_NP1728, fNAS_NP1728, 'r-', label = r'Np=1728, $\tilde{V}= 12^3$, mean(fNAS)=%4.2f (%%)'%mean_fNAS_NP1728)
plt.grid()
plt.axis([0, 50, 0, 12])
plt.xlabel(r'dimensionless time / $\tau_0$')
plt.ylabel('fraction of NAS (%)')
plt.legend(loc = 'lower right')
plt.show()

