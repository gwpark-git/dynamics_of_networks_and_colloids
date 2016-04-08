

from numpy import *
import matplotlib.pyplot as plt

mean_fNAS_dt01 = []

move_path='/Volumes/Task_REPO_MAC_PRO/works/Brownian_simulation/data/stochastic_simulation_repository/data/density_dependence'

fNAS0400_dt01 = 100.*loadtxt('%s/NP0400_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0400_LD10P3_C100.ener'%move_path)[:,4]/(400.*25.)
t_fNAS0400_dt01 = arange(size(fNAS0400_dt01))/100.
mean_fNAS0400_dt01 = mean(fNAS0400_dt01[size(fNAS0400_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0400_dt01)

fNAS0600_dt01 = 100.*loadtxt('%s/NP0600_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0600_LD10P3_C100.ener'%move_path)[:,4]/(600.*25.)
t_fNAS0600_dt01 = arange(size(fNAS0600_dt01))/100.
mean_fNAS0600_dt01 = mean(fNAS0600_dt01[size(fNAS0600_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0600_dt01)

fNAS0800_dt01 = 100.*loadtxt('%s/NP0800_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0800_LD10P3_C100.ener'%move_path)[:,4]/(800.*25.)
t_fNAS0800_dt01 = arange(size(fNAS0800_dt01))/100.
mean_fNAS0800_dt01 = mean(fNAS0800_dt01[size(fNAS0800_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0800_dt01)


fNAS1000_dt01 = 100.*loadtxt('%s/longer_computation/NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener'%move_path)[:,4]/(1000.*25.)
t_fNAS1000_dt01 = arange(size(fNAS1000_dt01))/100.
mean_fNAS1000_dt01 = mean(fNAS1000_dt01[size(fNAS1000_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1000_dt01)

fNAS1200_dt01 = 100.*loadtxt('%s/longer_computation/NP1200_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1200_LD10P3_C100.ener'%move_path)[:,4]/(1200.*25.)
t_fNAS1200_dt01 = arange(size(fNAS1200_dt01))/100.
mean_fNAS1200_dt01 = mean(fNAS1200_dt01[size(fNAS1200_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1200_dt01)

fNAS1400_dt01 = 100.*loadtxt('%s/longer_computation/NP1400_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1400_LD10P3_C100.ener'%move_path)[:,4]/(1400.*25.)
t_fNAS1400_dt01 = arange(size(fNAS1400_dt01))/100.
mean_fNAS1400_dt01 = mean(fNAS1400_dt01[size(fNAS1400_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1400_dt01)


# data before modify time step
# mean_fNDAS_dt01 = []
# fNDAS0400 = 100.*loadtxt('move_path/NP0400_LD10P3_RT100_block100/NP0400_LD10P3_C100.ener')[:,4]/(400.*25.)
# t_fNDAS0400 = arange(size(fNDAS0400))/1000.
# mean_fNDAS0400 = mean(fNDAS0400[size(fNDAS0400)/2:])
# mean_fNDAS.append(mean_fNDAS0400)

# fNDAS0600 = 100.*loadtxt('move_path/NP0600_LD10P3_RT100_NB100/NP0600_LD10P3_C100.ener')[:,4]/(600.*25.)
# t_fNDAS0600 = arange(size(fNDAS0600))/1000.
# mean_fNDAS0600 = mean(fNDAS0600[size(fNDAS0600)/2:])
# mean_fNDAS.append(mean_fNDAS0600)

# fNDAS0800 = 100.*loadtxt('move_path/NP0800_LD10P3_RT100_block100/NP0800_LD10P3_C100.ener')[:,4]/(800.*25.)
# t_fNDAS0800 = arange(size(fNDAS0800))/1000.
# mean_fNDAS0800 = mean(fNDAS0800[size(fNDAS0800)/2:])
# mean_fNDAS.append(mean_fNDAS0800)

# fNDAS1000 = 100.*loadtxt('../longer_computation/NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
# t_fNDAS1000 = arange(size(fNDAS1000))/1000.
# mean_fNDAS1000 = mean(fNDAS1000[size(fNDAS1000)/2:])
# mean_fNDAS.append(mean_fNDAS1000)

# fNDAS1200 = 100.*loadtxt('../longer_computation/NP1200_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1200_LD10P3_C100.ener')[:,4]/(1200.*25.)
# t_fNDAS1200 = arange(size(fNDAS1200))/1000.
# mean_fNDAS1200 = mean(fNDAS1200[size(fNDAS1200)/2:])
# mean_fNDAS.append(mean_fNDAS1200)

# fNDAS1000 = 100.*loadtxt('../longer_computation/NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
# t_fNDAS1000 = arange(size(fNDAS1000))/1000.
# mean_fNDAS1000 = mean(fNDAS1000[size(fNDAS1000)/2:])
# mean_fNDAS.append(mean_fNDAS1000)

# # data after modify time step
mean_fNDAS_dt01 = []

# move_path='/Volumes/Task_REPO_MAC_PRO/works/Brownian_simulation/data/stochastic_simulation/density_dependence'

fNDAS0400_dt01 = (100./2.)*loadtxt('NDAS_NP0400_RT100.dat')/(400.*25.)
for i,x in enumerate(fNDAS0400_dt01):
    if int(x) <= 0 and i > 1000:
        fNDAS0400_dt01 = fNDAS0400_dt01[:i]
t_fNDAS0400_dt01 = arange(size(fNDAS0400_dt01))/100.
mean_fNDAS0400_dt01 = mean(fNDAS0400_dt01[size(fNDAS0400_dt01)/2:])
mean_fNDAS_dt01.append(mean_fNDAS0400_dt01)

fNDAS0600_dt01 = (100./2.)*loadtxt('NDAS_NP0600_RT100.dat')/(600.*25.)
for i,x in enumerate(fNDAS0600_dt01):
    if int(x) <= 0 and i > 1000:
        fNDAS0600_dt01 = fNDAS0600_dt01[:i]
t_fNDAS0600_dt01 = arange(size(fNDAS0600_dt01))/100.
mean_fNDAS0600_dt01 = mean(fNDAS0600_dt01[size(fNDAS0600_dt01)/2:])
mean_fNDAS_dt01.append(mean_fNDAS0600_dt01)

fNDAS0800_dt01 = (100./2.)*loadtxt('NDAS_NP0800_RT100.dat')/(800.*25.)
for i,x in enumerate(fNDAS0800_dt01):
    if int(x) <= 0 and i > 1000:
        fNDAS0800_dt01 = fNDAS0800_dt01[:i]
t_fNDAS0800_dt01 = arange(size(fNDAS0800_dt01))/100.
mean_fNDAS0800_dt01 = mean(fNDAS0800_dt01[size(fNDAS0800_dt01)/2:])
mean_fNDAS_dt01.append(mean_fNDAS0800_dt01)

fNDAS1000_dt01 = (100./2.)*loadtxt('NDAS_NP1000_RT100.dat')/(1000.*25.)
for i,x in enumerate(fNDAS1000_dt01):
    if int(x) <= 0 and i > 1000:
        fNDAS1000_dt01 = fNDAS1000_dt01[:i]
t_fNDAS1000_dt01 = arange(size(fNDAS1000_dt01))/100.
mean_fNDAS1000_dt01 = mean(fNDAS1000_dt01[size(fNDAS1000_dt01)/2:])
mean_fNDAS_dt01.append(mean_fNDAS1000_dt01)

fNDAS1200_dt01 = (100./2.)*loadtxt('NDAS_NP1200_RT100.dat')/(1200.*25.)
for i,x in enumerate(fNDAS1200_dt01):
    if int(x) <= 0 and i > 1000:
        fNDAS1200_dt01 = fNDAS1200_dt01[:i]
t_fNDAS1200_dt01 = arange(size(fNDAS1200_dt01))/100.
mean_fNDAS1200_dt01 = mean(fNDAS1200_dt01[size(fNDAS1200_dt01)/2:])
mean_fNDAS_dt01.append(mean_fNDAS1200_dt01)

fNDAS1400_dt01 = (100./2.)*loadtxt('NDAS_NP1400_RT100.dat')/(1400.*25.)
for i,x in enumerate(fNDAS1400_dt01):
    if int(x) <= 0 and i > 1000:
        fNDAS1400_dt01 = fNDAS1400_dt01[:i]
t_fNDAS1400_dt01 = arange(size(fNDAS1400_dt01))/100.
mean_fNDAS1400_dt01 = mean(fNDAS1400_dt01[size(fNDAS1400_dt01)/2:])
mean_fNDAS_dt01.append(mean_fNDAS1400_dt01)

# fNDAS2000_dt01 = 100.*loadtxt('../longer_computation/NP2000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP2000_LD10P3_C100.ener')[:,4]/(2000.*25.)
# t_fNDAS2000_dt01 = arange(size(fNDAS2000_dt01))/100.
# mean_fNDAS2000_dt01 = mean(fNDAS2000_dt01[size(fNDAS2000_dt01)/2:])
# mean_fNDAS_dt01.append(mean_fNDAS2000_dt01)

# fNDAS4000_dt01 = 100.*loadtxt('../longer_computation/NP4000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP4000_LD10P3_C100.ener')[:,4]/(4000.*25.)
# t_fNDAS4000_dt01 = arange(size(fNDAS4000_dt01))/100.
# mean_fNDAS4000_dt01 = mean(fNDAS4000_dt01[size(fNDAS4000_dt01)/2:])
# mean_fNDAS_dt01.append(mean_fNDAS4000_dt01)



colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple', 'gray', 'black']
Np = [400, 600, 800, 1000, 1200, 1400]
NDP = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
# NVP = asarray(Np)*(4./3.)*pi*0.5**3/1000.
# plt.close()
# plt.figure(figsize=(6,8))
# plt.ion()
# plt.plot(t_fNDAS0400_dt01, fNDAS0400_dt01, '-', color=colP[0], label = 'NP=400, mean=%3.2f(%%)'%(mean_fNDAS0400_dt01))
# plt.plot(t_fNDAS0600_dt01, fNDAS0600_dt01, '-', color=colP[1], label = 'NP=600, mean=%3.2f(%%)'%(mean_fNDAS0600_dt01))
# plt.plot(t_fNDAS0800_dt01, fNDAS0800_dt01, '-', color=colP[2], label = 'NP=800, mean=%3.2f(%%)'%(mean_fNDAS0800_dt01))
# plt.plot(t_fNDAS1000_dt01, fNDAS1000_dt01, '-', color=colP[3], label = 'NP=1000, mean=%3.2f(%%)'%(mean_fNDAS1000_dt01))
# plt.plot(t_fNDAS1200_dt01, fNDAS1200_dt01, '-', color=colP[4], label = 'NP=1200, mean=%3.2f(%%)'%(mean_fNDAS1200_dt01))
# plt.plot(t_fNDAS1400_dt01, fNDAS1400_dt01, '-', color=colP[5], label = 'NP=1400, mean=%3.2f(%%)'%(mean_fNDAS1400_dt01))

# plt.grid()
# plt.yticks(mean_fNDAS_dt01, ['%3.2f' % val for val in mean_fNDAS_dt01])
# # plt.xticks(Np)
# plt.axis([0, 50, 0, 25])

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
# fontP.set_size('small')
# plt.legend(loc = 'lower right', ncol=1, prop=fontP)
# plt.xlabel('topological dimensionless time')
# plt.ylabel('fraction of NDAS (%)')
# plt.show()


NP = [400, 600, 800, 1000, 1200, 1400]
# NP = [10, 15, 20, 25, 30, 35]
# NP_dat=NP[:-2]
# mean_fNDAS_dt01_dat = mean_fNDAS_dt01[:-2]
# mean_fNDAS_dt01 = asarange([])
plt.close()
plt.figure(figsize=(8,6))
plt.ion()
ax = plt.subplot(111)
ax.plot(NDP, mean_fNAS_dt01, 'bo-', markerfacecolor='black', label = 'including multiple connections')
ax.plot(NDP, mean_fNDAS_dt01, 'rs-', markerfacecolor='white', label = 'excluding multiple connections')
text_legend_fNAS = []
text_legend_fNDAS = []
for i in range(size(NP)):
    # ax.plot(NP[i], mean_fNAS_dt01[i], 'o', markersize=6, markerfacecolor=colP[i], label = 'NP=%d'%(NP[i]))
    ax.plot(NDP[i], mean_fNAS_dt01[i], 'o', markersize=8, markerfacecolor=colP[i])
    text_legend_fNAS.append('%3.1f'%(mean_fNAS_dt01[i]))
    ax.annotate('w=%3.2f'%(mean_fNAS_dt01[i]/mean_fNDAS_dt01[i]), xy=(NDP[i]+0.05, mean_fNAS_dt01[i]))
    ax.plot([0.3, NDP[i]], [mean_fNAS_dt01[i], mean_fNAS_dt01[i]], 'b:', alpha=0.5)
for i in range(size(NP)):
    # ax.plot(NP[i], mean_fNDAS_dt01[i], 'o', markersize=6, markerfacecolor=colP[i], label = 'NP=%d'%(NP[i]))
    ax.plot(NDP[i], mean_fNDAS_dt01[i], 's', markersize=8, markeredgecolor=colP[i], markerfacecolor='white')
    # text_legend_fNDAS.append('%3.1f\nf=%3.1f'%(mean_fNDAS_dt01[i], (mean_fNDAS_dt01[i]/100.)*25.*2.))
    text_legend_fNDAS.append('%3.1f'%(mean_fNDAS_dt01[i]))
    # ax.annotate('f=%3.1f'%((mean_fNDAS_dt01[i]/100.)*25.*2.), xy=(NDP[i]+0.05, mean_fNDAS_dt01[i]))
    ax.annotate('f=%3.1f'%((mean_fNDAS_dt01[i]/100.)*25.*2.), xy=(NDP[i]+0.05, mean_fNDAS_dt01[i]))
    ax.plot([NDP[i], 1.6], [mean_fNDAS_dt01[i], mean_fNDAS_dt01[i]], 'r:', alpha=0.5)
    
ax.legend(loc = 'upper left', numpoints=1, prop=fontP)
ax.set_xlabel(r'micelle number density, $\nu_m$')
ax.set_ylabel('bridge chain fraction (%)')
ax.set_xticks(NDP)
ax.set_yticks(mean_fNAS_dt01)
ax.set_yticklabels(text_legend_fNAS)

# ax.grid(axis='x')
# ax.grid(color='blue', axis='y', fillstyle='left')
ax.axis([0.3, 1.6, 2, 20])
ax2 = ax.twinx()
ax2.axis([0.3, 1.6, 2, 20])
ax2.set_yticks(mean_fNDAS_dt01)
ax2.set_yticklabels(text_legend_fNDAS)
# ax2.set_ylabel('bridge chain fraction(%)')
# ax2.grid(color='red', axis='y', fillstyle='right')
# ax2.axis([300, 1500, 
plt.show()
