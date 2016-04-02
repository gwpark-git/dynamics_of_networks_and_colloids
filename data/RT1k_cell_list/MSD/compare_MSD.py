
from numpy import *
import matplotlib.pyplot as plt
from scipy.stats import linregress

MSD_EQ_NP0400 = loadtxt('MSD_EQ_NP0400_half.dat')[:,1]
Nt_EQ_NP0400 = size(MSD_EQ_NP0400)
t_EQ_NP0400 = arange(Nt_EQ_NP0400)/1000.
reg_EQ_NP0400 = linregress(t_EQ_NP0400[Nt_EQ_NP0400/2:], MSD_EQ_NP0400[Nt_EQ_NP0400/2:])

MSD_EQ_NP1000 = loadtxt('MSD_EQ_NP1000_half.dat')[:,1]
Nt_EQ_NP1000 = size(MSD_EQ_NP1000)
t_EQ_NP1000 = arange(Nt_EQ_NP1000)/1000.
reg_EQ_NP1000 = linregress(t_EQ_NP1000[Nt_EQ_NP1000/2:], MSD_EQ_NP1000[Nt_EQ_NP1000/2:])



MSD_NP0400 = loadtxt('MSD_NP0400.dat')[:,1]
Nt_NP0400 = size(MSD_NP0400)
t_NP0400 = arange(Nt_NP0400)/100.
reg_NP0400 = linregress(t_NP0400[Nt_NP0400/2:], MSD_NP0400[Nt_NP0400/2:])


MSD_NP1000 = loadtxt('MSD_NP1000.dat')[:,1]
Nt_NP1000 = size(MSD_NP1000)
t_NP1000 = arange(Nt_NP1000)/100.
reg_NP1000 = linregress(t_NP1000[Nt_NP1000/2:], MSD_NP1000[Nt_NP1000/2:])


t_theory = linspace(0, 100, 100)
y_theory = 6*t_theory*10
plt.close()
plt.ion()
plt.plot(t_theory, y_theory, 'k--', label = r'Theory (pure Brownian), slope=$6Rt=60$, $\tilde{D} = 10 \Rightarrow t_B = 1/\tilde{D} = 0.10$')
plt.plot(t_EQ_NP0400, MSD_EQ_NP0400, 'b-', linewidth=2, alpha=0.5, label = r'Np=0400, slope=%4.1f, $\tilde{D}$=%4.2f, $t_B$=%4.2f'%(reg_EQ_NP0400[0], reg_EQ_NP0400[0]/6., 1/(reg_EQ_NP0400[0]/6.)))
plt.plot(t_EQ_NP1000, MSD_EQ_NP1000, 'r-', linewidth=2, alpha=0.5, label = r'Np=0400, slope=%4.1f, $\tilde{D}$=%4.2f, $t_B$=%4.2f'%(reg_EQ_NP1000[0], reg_EQ_NP1000[0]/6., 1/(reg_EQ_NP1000[0]/6.)))

plt.plot(t_NP0400, MSD_NP0400, 'b-', label = r'Np=0400, slope=%4.1f, $\tilde{D}$=%4.2f, $t_B$=%4.2f'%(reg_NP0400[0], reg_NP0400[0]/6., 1/(reg_NP0400[0]/6.)))
plt.plot(t_NP1000, MSD_NP1000, 'r-', label = r'Np=1000, slope=%4.1f, $\tilde{D}$=%4.2f, $t_B$=%4.2f'%(reg_NP1000[0], reg_NP1000[0]/6., 1/(reg_NP1000[0]/6.)))

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'upper left', prop=fontP)
plt.axis([0, 0.5, 0, 20])

plt.xlabel(r'dimensionless time / $\tau_0$')
plt.ylabel(r'MSD, $\langle (\tilde{\mathbf{p}}_i(t) - \tilde{\mathbf{p}}_i(0))^2\rangle_i$')

plt.grid()
plt.show()
