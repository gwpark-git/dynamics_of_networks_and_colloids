
from numpy import *
import matplotlib.pyplot as plt
from scipy.stats import linregress

MSD_EQ_NP0800 = loadtxt('MSD_NP0800_EQ.dat')[:,1]
Nt_EQ_NP0800 = size(MSD_EQ_NP0800)
t_EQ_NP0800 = arange(Nt_EQ_NP0800)/100.
reg_EQ_NP0800 = linregress(t_EQ_NP0800[Nt_EQ_NP0800/2:], MSD_EQ_NP0800[Nt_EQ_NP0800/2:])


# MSD_NP0400_RT1 = loadtxt('MSD_NP0400_RT1.dat')[:,1]
# Nt_NP0400_RT1 = size(MSD_NP0400_RT1)
# t_NP0400_RT1 = arange(Nt_NP0400_RT1)/1000.
# reg_NP0400_RT1 = linregress(t_NP0400_RT1[Nt_NP0400_RT1/2:], MSD_NP0400_RT1[Nt_NP0400_RT1/2:])


# MSD_NP1000_RT1 = loadtxt('MSD_NP1000_RT1.dat')[:,1]
# Nt_NP1000_RT1 = size(MSD_NP1000_RT1)
# t_NP1000_RT1 = arange(Nt_NP1000_RT1)/1000.
# reg_NP1000_RT1 = linregress(t_NP1000_RT1[Nt_NP1000_RT1/2:], MSD_NP1000_RT1[Nt_NP1000_RT1/2:])


t_theory = linspace(0, 100, 100)
y_theory = 6*t_theory*10
plt.close()
plt.ion()
plt.plot(t_theory, y_theory, 'k:', linewidth=3, alpha=0.3, label = r'Theory for pure Brownian: slope=$6R_t=60$, $\tilde{D} = 10 \Rightarrow \tau_c = 1/\tilde{D} = 0.10 \tau_0$')
plt.plot(t_EQ_NP0800, MSD_EQ_NP0800, 'b-', linewidth=3, alpha=0.3, label = r'Np=0400 (repulsive Brownian), slope=%4.1f, $\tilde{D}$=%4.2f, $\tau_c$=%4.2f $\tau_0$'%(reg_EQ_NP0800[0], reg_EQ_NP0800[0]/6., 1/(reg_EQ_NP0800[0]/6.)))


# plt.plot(t_NP0400_RT1, MSD_NP0400_RT1, 'b--', label = r'Np=0400 (Rt=1), slope=%4.1f, $\tilde{D}$=%4.2f, $\tau_c$=%4.2f $\tau_0$'%(reg_NP0400_RT1[0], reg_NP0400_RT1[0]/6., 1/(reg_NP0400_RT1[0]/6.)))
# plt.plot(t_NP1000_RT1, MSD_NP1000_RT1, 'r--', label = r'Np=1000 (Rt=1), slope=%4.1f, $\tilde{D}$=%4.2f, $\tau_c$=%4.2f $\tau_0$'%(reg_NP1000_RT1[0], reg_NP1000_RT1[0]/6., 1/(reg_NP1000_RT1[0]/6.)))

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'upper left', prop=fontP)
plt.axis([0, 10, 0, 400])

plt.xlabel(r'dimensionless time / $\tau_0$')
plt.ylabel(r'MSD, $\langle (\tilde{\mathbf{p}}_i(t) - \tilde{\mathbf{p}}_i(0))^2\rangle_i$')

plt.grid()
plt.show()
