
# from numpy import *
# import matplotlib.pyplot as plt

# import sys
# sys.path.append('../../post_processing')
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
t_peak_NP0400, acf_peak_NP0400 = t_acf_NP0400[90:100][argmax(acf_NP0400[90:100])], max(acf_NP0400[90:100])
peak_info.append([10, t_peak_NP0400, acf_peak_NP0400, acf_NP0400[0]])
t_peak_NP0600, acf_peak_NP0600 = t_acf_NP0600[100:120][argmax(acf_NP0600[100:120])], max(acf_NP0600[100:120])
peak_info.append([15, t_peak_NP0600, acf_peak_NP0600, acf_NP0600[0]])
t_peak_NP0800, acf_peak_NP0800 = t_acf_NP0800[100:160][argmax(acf_NP0800[100:160])], max(acf_NP0800[100:160])
peak_info.append([20, t_peak_NP0800, acf_peak_NP0800, acf_NP0800[0]])
t_peak_NP1000, acf_peak_NP1000 = t_acf_NP1000[200:430][argmax(acf_NP1000[200:430])], max(acf_NP1000[200:430])
peak_info.append([25, t_peak_NP1000, acf_peak_NP1000, acf_NP1000[0]])
t_peak_NP1200, acf_peak_NP1200 = t_acf_NP1200[600:700][argmax(acf_NP1200[600:700])], max(acf_NP1200[600:700])
peak_info.append([30, t_peak_NP1200, acf_peak_NP1200, acf_NP1200[0]])
t_peak_NP1400, acf_peak_NP1400 = t_acf_NP1400[1200:1400][argmax(acf_NP1400[1200:1400])], max(acf_NP1400[1200:1400])
peak_info.append([35, t_peak_NP1400, acf_peak_NP1400, acf_NP1400[0]])
peak_info = asarray(peak_info)

plt.close()
plt.ion()
plt.figure(figsize=(6,8))
ax = plt.subplot(111)
ax2 = ax.twinx()
leg_1 = []
leg_2 = []
tmp_leg_1, = ax.loglog(peak_info[:,0], peak_info[:,1], 'b-', label = r'$\tau_{eff}/\tau_0$ (left axis)')
leg_1.append(tmp_leg_1)
tmp_leg_2, = ax2.loglog(peak_info[:,0], peak_info[:,3], 'r-', label = r'$\tilde{G}(0)$ (right axis)')
leg_1.append(tmp_leg_2)

nuR3 = [10, 15, 20, 25, 30, 35]
f = [1.5, 2.1, 2.9, 3.7, 4.6, 5.6]
w = [1.54, 1.55, 1.56, 1.58, 1.60, 1.64]
cnt=0
tmp_leg_1, = ax.loglog(peak_info[cnt,0], peak_info[cnt,1],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt], label = r'$\nu_0=%d/R_0^3$, f=%4.1f, w=%4.2f'%(nuR3[cnt], f[cnt], w[cnt]))
leg_1.append(tmp_leg_1)

tmp_leg_2, = ax2.loglog(peak_info[cnt,0], peak_info[cnt,3], linestyle='none',  marker=symP[cnt], markerfacecolor=colP[cnt])
leg_2.append(tmp_leg_2)

cnt+=1
tmp_leg_1, = ax.loglog(peak_info[cnt,0], peak_info[cnt,1],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt], label = r'$\nu_0=%d/R_0^3$, f=%4.1f, w=%4.2f'%(nuR3[cnt], f[cnt], w[cnt]))
leg_1.append(tmp_leg_1)

tmp_leg_2, = ax2.loglog(peak_info[cnt,0], peak_info[cnt,3],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt])
leg_2.append(tmp_leg_2)

cnt+=1 
tmp_leg_1, = ax.loglog(peak_info[cnt,0], peak_info[cnt,1],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt], label = r'$\nu_0=%d/R_0^3$, f=%4.1f, w=%4.2f'%(nuR3[cnt], f[cnt], w[cnt]))
leg_1.append(tmp_leg_1)
tmp_leg_2, = ax2.loglog(peak_info[cnt,0], peak_info[cnt,3],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt])
leg_2.append(tmp_leg_2)
cnt+=1 
tmp_leg_1, = ax.loglog(peak_info[cnt,0], peak_info[cnt,1],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt], label = r'$\nu_0=%d/R_0^3$, f=%4.1f, w=%4.2f'%(nuR3[cnt], f[cnt], w[cnt]))
leg_1.append(tmp_leg_1)
tmp_leg_2, = ax2.loglog(peak_info[cnt,0], peak_info[cnt,3],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt])
leg_2.append(tmp_leg_2)

cnt+=1 
tmp_leg_1, = ax.loglog(peak_info[cnt,0], peak_info[cnt,1],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt], label = r'$\nu_0=%d/R_0^3$, f=%4.1f, w=%4.2f'%(nuR3[cnt], f[cnt], w[cnt]))
leg_1.append(tmp_leg_1)
tmp_leg_2, = ax2.loglog(peak_info[cnt,0], peak_info[cnt,3],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt])
leg_2.append(tmp_leg_2)

cnt+=1 
tmp_leg_1, = ax.loglog(peak_info[cnt,0], peak_info[cnt,1],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt], label = r'$\nu_0=%d/R_0^3$, f=%4.1f, w=%4.2f'%(nuR3[cnt], f[cnt], w[cnt]))
leg_1.append(tmp_leg_1)
tmp_leg_2, = ax2.loglog(peak_info[cnt,0], peak_info[cnt,3],  linestyle='none', marker=symP[cnt], markerfacecolor=colP[cnt])
leg_2.append(tmp_leg_2)

ax.axis([9, 36, 0.8, 16])
ax.set_xticks(peak_info[:,0])
ax.set_xticklabels(['%4.0f'% val for val in peak_info[:,0]])
ax.set_xlabel(r'number density of chains, $\nu_0R_0^3$')
ax.set_ylabel(r'Effective network relaxation time / loop relaxation time, $\tau_{eff}/\tau_0$')
ax2.set_ylabel(r'$\tilde{G}(0)$')
ax.grid(axis='x')
# leg_list = [legs for sublist in [leg_1, leg_2] for legs in sublist]
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
ax.legend(loc = 'upper left', handles=leg_1, prop=fontP, numpoints=1, ncol=1)
# ax3 = ax.twiny()
# ax3.loglog(100, 100, 'wo')
# ax3.set_xticks(peak_info[:,0])
# ax3.set_xticklabels(['f=%4.1f'% val for val in [1.5, 2.1, 2.9, 3.7, 4.6, 5.6]])
# ax3.axis([0, 36, 0.8, 16])
# ax3.set_xlabel('f')

# plt.close()
# plt.ion()
# plt.figure(figsize=(11,6))
# plt.loglog(t_acf_NP0400[:100], acf_NP0400[:100], '-', color=colP[0], label = r'$\nu_0=10/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP0400[0])
# plt.loglog(t_acf_NP0600[:120], acf_NP0600[:120], '-', color=colP[1], label = r'$\nu_0=15/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP0600[0])

# plt.loglog(t_acf_NP0800[:160], acf_NP0800[:160], '-', color=colP[2], label = r'$\nu_0=20/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP0800[0])

# plt.loglog(t_acf_NP1000[:430], acf_NP1000[:430], '-', color=colP[3], label = r'$\nu_0=25/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP1000[0])

# plt.loglog(t_acf_NP1200[:700], acf_NP1200[:700], '-', color=colP[4], label = r'$\nu_0=30/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP1200[0])

# plt.loglog(t_acf_NP1400[:1400], acf_NP1400[:1400], '-', color=colP[5], label = r'$\nu_0=35/R_0^3$, $\tilde{G}(0)$=%4.1e'%acf_NP1400[0])

# plt.loglog(t_peak_NP0400, acf_peak_NP0400, marker=symP[0], color=colP[0], label = r'$\tau_{eff}$=%4.2f $\tau_0$, $\tilde{G}(\tau_{eff})$=%4.1e'%(t_peak_NP0400, acf_peak_NP0400))
# plt.loglog(t_peak_NP0600, acf_peak_NP0600, marker=symP[1], color=colP[1], label = r'$\tau_{eff}$=%4.2f $\tau_0$, $\tilde{G}(\tau_{eff})$=%4.1e'%(t_peak_NP0600, acf_peak_NP0600))
# plt.loglog(t_peak_NP0800, acf_peak_NP0800, marker=symP[2], color=colP[2], label = r'$\tau_{eff}$=%4.2f $\tau_0$, $\tilde{G}(\tau_{eff})$=%4.1e'%(t_peak_NP0800, acf_peak_NP0800))
# plt.loglog(t_peak_NP1000, acf_peak_NP1000, marker=symP[3], color=colP[3], label = r'$\tau_{eff}$=%4.2f $\tau_0$, $\tilde{G}(\tau_{eff})$=%4.1e'%(t_peak_NP1000, acf_peak_NP1000))
# plt.loglog(t_peak_NP1200, acf_peak_NP1200, marker=symP[4], color=colP[4], label = r'$\tau_{eff}$=%4.2f $\tau_0$, $\tilde{G}(\tau_{eff})$=%4.1e'%(t_peak_NP1200, acf_peak_NP1200))
# plt.loglog(t_peak_NP1400, acf_peak_NP1400, marker=symP[5], color=colP[5], label = r'$\tau_{eff}$=%4.2f $\tau_0$, $\tilde{G}(\tau_{eff})$=%4.1e'%(t_peak_NP1400, acf_peak_NP1400))


# plt.xlabel(r'topological time (normalized by $\tau_0$)')
# plt.ylabel('dimensionless stress autocorrelation for bridges')
# plt.axis([0.1, 20, 10**-8, 10**-2])

# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# fontP.set_size('x-small')
# plt.grid()

# plt.legend(loc = 'lower left', numpoints=1, ncol=2, prop=fontP)
# # plt.plot(t_sf_NP1000, sf_NP1000, 'b-')

# plt.axes([.55, .2, .3, .3])
# # peak_info[:,0] = peak_info[:,0]**(2./3.)
# plt.loglog(peak_info[:,0], peak_info[:,1], 'k-')
# cnt=0
# plt.loglog(peak_info[cnt,0], peak_info[cnt,1], marker=symP[cnt], color=colP[cnt])
# cnt+=1
# plt.loglog(peak_info[cnt,0], peak_info[cnt,1], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(peak_info[cnt,0], peak_info[cnt,1], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(peak_info[cnt,0], peak_info[cnt,1], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(peak_info[cnt,0], peak_info[cnt,1], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(peak_info[cnt,0], peak_info[cnt,1], marker=symP[cnt], color=colP[cnt])
# plt.axis([9, 36, 0.8, 16])
# plt.xticks(peak_info[:,0], ['%4.0f'% val for val in peak_info[:,0]])
# plt.xlabel(r'number density of chains, $\nu_0R_0^3$')
# plt.ylabel(r'          $\tau_{eff}/\tau_0$')
# plt.grid(axis='x')
# # plt.loglog(
# plt.show()
