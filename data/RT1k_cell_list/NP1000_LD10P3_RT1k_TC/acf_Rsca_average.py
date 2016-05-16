
from numpy import *
import matplotlib.pyplot as plt

Np      = 1000
Nc      = 10

N_tot   = Nc*Np

t = arange(50*100)/100.
acf_av_Rsca = zeros(size(t))
for k in range(N_tot):
    acf_Rsca_k = loadtxt('tracking_individual_chain/acf_Rsca_%06d.dat'%k)[:,1]
    acf_av_Rsca[:size(acf_Rsca_k)] += acf_Rsca_k
acf_av_Rsca /= float(N_tot)

ref_fraction = asarray([[-1, 0.1041],
                        [21, 0.1041]])

from scipy.stats import linregress
val_reg = linregress(t[:100], log(acf_av_Rsca[:100]))
y_reg = exp(val_reg[0]*t)

plt.close()
plt.ion()
plt.figure(figsize=(6, 6))
plt.plot(t, acf_av_Rsca/acf_av_Rsca[0], 'r-', linewidth=2, label = 'Np=1000')
plt.plot(t, y_reg, 'r--', label = r'reg up to 1$\tau_0$, tc=%3.2f$\tau_0$'%(-1./val_reg[0]))
plt.plot(ref_fraction[:,0], ref_fraction[:,1], 'k:', linewidth=3, label = 'fNAS=10.41%')
plt.fill_between(t, acf_av_Rsca/acf_av_Rsca[0], 0, color='red', alpha=0.2)
plt.legend(loc = 'upper right')
plt.grid()
plt.xlabel(r'dimensionless time, $\beta_0 t$')
plt.ylabel('averaged autocorrelation')
plt.axis([-1, 21, -0.05, 1.2])
plt.show()
