

from numpy import *
import matplotlib.pyplot as plt

# data before modify time step
mean_fNAS_dt01 = []

fNAS0400 = 100.*loadtxt('../NP0400_LD10P3_RT100_block100/NP0400_LD10P3_C100.ener')[:,4]/(400.*25.)
t_fNAS0400 = arange(size(fNAS0400))/1000.
mean_fNAS0400 = mean(fNAS0400[size(fNAS0400)/2:])
mean_fNAS.append(mean_fNAS0400)

fNAS0600 = 100.*loadtxt('../NP0600_LD10P3_RT100_NB100/NP0600_LD10P3_C100.ener')[:,4]/(600.*25.)
t_fNAS0600 = arange(size(fNAS0600))/1000.
mean_fNAS0600 = mean(fNAS0600[size(fNAS0600)/2:])
mean_fNAS.append(mean_fNAS0600)

fNAS0800 = 100.*loadtxt('../NP0800_LD10P3_RT100_block100/NP0800_LD10P3_C100.ener')[:,4]/(800.*25.)
t_fNAS0800 = arange(size(fNAS0800))/1000.
mean_fNAS0800 = mean(fNAS0800[size(fNAS0800)/2:])
mean_fNAS.append(mean_fNAS0800)

fNAS1000 = 100.*loadtxt('../NP1000_LD10P3_RT100_NB100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
t_fNAS1000 = arange(size(fNAS1000))/1000.
mean_fNAS1000 = mean(fNAS1000[size(fNAS1000)/2:])
mean_fNAS.append(mean_fNAS1000)

# data after modify time step
mean_fNAS_dt01 = []

fNAS0400_dt01 = 100.*loadtxt('../NP0400_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0400_LD10P3_C100.ener')[:,4]/(400.*25.)
t_fNAS0400_dt01 = arange(size(fNAS0400_dt01))/100.
mean_fNAS0400_dt01 = mean(fNAS0400_dt01[size(fNAS0400_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0400_dt01)

fNAS0600_dt01 = 100.*loadtxt('../NP0600_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0600_LD10P3_C100.ener')[:,4]/(600.*25.)
t_fNAS0600_dt01 = arange(size(fNAS0600_dt01))/100.
mean_fNAS0600_dt01 = mean(fNAS0600_dt01[size(fNAS0600_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0600_dt01)

fNAS0800_dt01 = 100.*loadtxt('../NP0800_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0800_LD10P3_C100.ener')[:,4]/(800.*25.)
t_fNAS0800_dt01 = arange(size(fNAS0800_dt01))/100.
mean_fNAS0800_dt01 = mean(fNAS0800_dt01[size(fNAS0800_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0800_dt01)

fNAS1000_dt01 = 100.*loadtxt('../NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
t_fNAS1000_dt01 = arange(size(fNAS1000_dt01))/100.
mean_fNAS1000_dt01 = mean(fNAS1000_dt01[size(fNAS1000_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1000_dt01)


colP = ['blue', 'green', 'brown', 'red']
Np = [400, 500, 600, 700]
plt.close()
plt.figure(figsize=(11,6))
plt.ion()
# plot data before time step modification
plt.plot(t_fNAS0400, fNAS0400, '-', color=colP[0], linewidth=3, alpha=0.2, label = 'Np=400, dt=1e-3, mean=%3.2f(%%)'%(mean_fNAS0400))
plt.plot(t_fNAS0600, fNAS0600, '-', color=colP[1], linewidth=3, alpha=0.2, label = 'Np=600, dt=1e-3, mean=%3.2f(%%)'%(mean_fNAS0600))
plt.plot(t_fNAS0800, fNAS0800, '-', color=colP[2], linewidth=3, alpha=0.2, label = 'Np=800, dt=1e-3, mean=%3.2f(%%)'%(mean_fNAS0800))
plt.plot(t_fNAS1000, fNAS1000, '-', color=colP[3], linewidth=3, alpha=0.2, label = 'Np=1000, dt=1e-3, last=%3.2f(%%)'%(fNAS1000[-1]))

# plot data after time step modification
plt.plot(t_fNAS0400_dt01, fNAS0400_dt01, '-', color=colP[0], label = 'dt=1e-2, mean=%3.2f(%%), diff=%3.2f(%%)'%(mean_fNAS0400_dt01, abs(100.*mean_fNAS0400/mean_fNAS0400_dt01 -100.)))
plt.plot(t_fNAS0600_dt01, fNAS0600_dt01, '-', color=colP[1], label = 'dt=1e-2, mean=%3.2f(%%), diff=%3.2f(%%)'%(mean_fNAS0600_dt01, abs(100.*mean_fNAS0600/mean_fNAS0600_dt01 - 100.)))
plt.plot(t_fNAS0800_dt01, fNAS0800_dt01, '-', color=colP[2], label = 'dt=1e-2, mean=%3.2f(%%), diff=%3.2f(%%)'%(mean_fNAS0800_dt01, abs(100.*mean_fNAS0800/mean_fNAS0800_dt01 - 100.)))
plt.plot(t_fNAS1000_dt01, fNAS1000_dt01, '-', color=colP[3], label = 'dt=1e-2, mean=%3.2f(%%), diff=%3.2f(%%)'%(mean_fNAS1000_dt01, abs(100.*fNAS1000[-1]/mean_fNAS1000_dt01 - 100)))

plt.grid()
plt.yticks(mean_fNAS_dt01, ['%3.2f' % val for val in mean_fNAS_dt01])
# plt.xticks(Np)
plt.axis([0, 50, 0, 13])

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'lower right', ncol=2, prop=fontP)
plt.xlabel('topological dimensionless time')
plt.ylabel('fraction of NAS (%)')
plt.show()

