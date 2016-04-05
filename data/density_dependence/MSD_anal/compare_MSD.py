
from numpy import *
import matplotlib.pyplot as plt
from scipy.stats import linregress


index_to_tau_B_RT1 = 100.
index_to_tau_0_RT1 = 100.


MSD_NP0400 = loadtxt('MSD_NP0400.dat')[:,1]
Nt_NP0400 = size(MSD_NP0400)
t_NP0400 = arange(Nt_NP0400)/index_to_tau_B_RT1
reg_NP0400 = linregress(t_NP0400[Nt_NP0400/2:], MSD_NP0400[Nt_NP0400/2:])

MSD_NP0600 = loadtxt('MSD_NP0600.dat')[:,1]
Nt_NP0600 = size(MSD_NP0600)
t_NP0600 = arange(Nt_NP0600)/index_to_tau_B_RT1
reg_NP0600 = linregress(t_NP0600[Nt_NP0600/2:], MSD_NP0600[Nt_NP0600/2:])

MSD_NP0800 = loadtxt('MSD_NP0800.dat')[:,1]
Nt_NP0800 = size(MSD_NP0800)
t_NP0800 = arange(Nt_NP0800)/index_to_tau_B_RT1
reg_NP0800 = linregress(t_NP0800[Nt_NP0800/2:], MSD_NP0800[Nt_NP0800/2:])



MSD_NP1000 = loadtxt('MSD_NP1000.dat')[:,1]
Nt_NP1000 = size(MSD_NP1000)
t_NP1000 = arange(Nt_NP1000)/index_to_tau_B_RT1
reg_NP1000 = linregress(t_NP1000[Nt_NP1000/2:], MSD_NP1000[Nt_NP1000/2:])

MSD_NP1200 = loadtxt('MSD_NP1200.dat')[:,1]
Nt_NP1200 = size(MSD_NP1200)
t_NP1200 = arange(Nt_NP1200)/index_to_tau_B_RT1
reg_NP1200 = linregress(t_NP1200[Nt_NP1200/2:], MSD_NP1200[Nt_NP1200/2:])

MSD_NP1400 = loadtxt('MSD_NP1400.dat')[:,1]
Nt_NP1400 = size(MSD_NP1400)
t_NP1400 = arange(Nt_NP1400)/index_to_tau_B_RT1
reg_NP1400 = linregress(t_NP1400[Nt_NP1400/2:], MSD_NP1400[Nt_NP1400/2:])



colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple']
t_theory = linspace(0, 100, 100)
y_theory = 6*t_theory
plt.close()
plt.figure(figsize=(11,6))
plt.ion()
plt.plot(t_theory, y_theory, 'k:', linewidth=3, alpha=0.3, label = r'Theory for pure Brownian: slope=$6$, $\tilde{D} = 1 \Rightarrow \tau_c$ = 1.00 $\tau_B$')

cnt = 0
plt.plot(t_NP0400, MSD_NP0400, '-', color=colP[cnt], label = r'Np=0400 (Rt=1, Nc=25, fNAS=4.6%%), slope=%4.1f, $\tilde{D}$=%4.3f, $\tau_c$=%4.3f $\tau_B$'%(reg_NP0400[0], reg_NP0400[0]/6., 1/(reg_NP0400[0]/6.)))
cnt += 1
plt.plot(t_NP0600, MSD_NP0600, '-', color=colP[cnt], label = r'Np=0600 (Rt=1, Nc=25, fNAS=6.6%%), slope=%4.1f, $\tilde{D}$=%4.3f, $\tau_c$=%4.3f $\tau_B$'%(reg_NP0600[0], reg_NP0600[0]/6., 1/(reg_NP0600[0]/6.)))
cnt += 1
plt.plot(t_NP0800, MSD_NP0800, '-', color=colP[cnt], label = r'Np=0800 (Rt=1, Nc=25, fNAS=8.9%%), slope=%4.1f, $\tilde{D}$=%4.3f, $\tau_c$=%4.3f $\tau_B$'%(reg_NP0800[0], reg_NP0800[0]/6., 1/(reg_NP0800[0]/6.)))
cnt += 1
plt.plot(t_NP1000, MSD_NP1000, '-', color=colP[cnt], label = r'Np=1000 (Rt=1, Nc=25, fNAS=11.6%%), slope=%4.1f, $\tilde{D}$=%4.3f, $\tau_c$=%4.3f $\tau_B$'%(reg_NP1000[0], reg_NP1000[0]/6., 1/(reg_NP1000[0]/6.)))
cnt += 1
plt.plot(t_NP1200, MSD_NP1200, '-', color=colP[cnt], label = r'Np=1200 (Rt=1, Nc=25, fNAS=14.8%%), slope=%4.1f, $\tilde{D}$=%4.3f, $\tau_c$=%4.3f $\tau_B$'%(reg_NP1200[0], reg_NP1200[0]/6., 1/(reg_NP1200[0]/6.)))
cnt += 1
plt.plot(t_NP1400, MSD_NP1400, '-', color=colP[cnt], label = r'Np=1400 (Rt=1, Nc=25, fNAS=18.3%%), slope=%4.1f, $\tilde{D}$=%4.3f, $\tau_c$=%4.3f $\tau_B$'%(reg_NP1400[0], reg_NP1400[0]/6., 1/(reg_NP1400[0]/6.)))


from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'upper left', prop=fontP)
plt.xticks(range(5))
plt.yticks(range(0, 31, 6))

plt.axis([0, 5, 0, 12])


plt.xlabel(r'dimensionless time / $\tau_B$')
plt.ylabel(r'MSD, $\langle (\tilde{\mathbf{p}}_i(t) - \tilde{\mathbf{p}}_i(0))^2\rangle_i$')

plt.grid()
plt.show()

# plt.close()
# plt.ion()
# plt.plot(t_NP1200, MSD_NP1200, 'b-')
# plt.grid()
# plt.show()
