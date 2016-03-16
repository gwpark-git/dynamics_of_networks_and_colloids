
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
fp = FontProperties()
fp.set_size('x-small')

file_list = ['NP1350_LD15P3_C100_SF15_RT010.NAS',
             'NP1350_LD15P3_C100_SF15_RT050.NAS',
             'NP1350_LD15P3_C100_SF15_RT100.NAS',
             'NP1350_LD15P3_C100_SF15_RT200.NAS',
             'NP1350_LD15P3_C100_SF15_RT300.NAS',
             'NP1350_LD15P3_C100_SF15_RT500.NAS',
             'NP1350_LD15P3_C100_SF15_RT1k.NAS']
dt = 0.001
RT = asarray([10., 50., 100., 200., 300., 500., 1000.])
delta_t = RT*dt
Nf = size(file_list)
dat = []
for fn in file_list:
    dat.append(loadtxt(fn))
Ntot = 1350*25.
av = zeros(Nf)
for i in range(Nf):
    dat[i][:, 0] -= dat[i][0, 0]
    Ni = shape(dat[i])[0]
    av[i] = mean(dat[i][Ni/2:,1])

csel = ['blue', 'cyan', 'red', 'purple', 'green', 'brown', 'black']



plt.clf()
plt.ion()
plt.figure(figsize=(11, 6))
ref_line = asarray([[0, mean(av[:3])],
                    [1, mean(av[:3])]])
ref_cutoff = asarray([[0.2, 1500],
                      [0.2, 2000]])
plt.plot(delta_t, av, 'bo-', markerfacecolor='white', markeredgecolor='blue', label = 'NAS vs topological time step')
plt.plot(ref_line[:,0], ref_line[:,1], 'k-', label= 'ref line, %4.2f(%6.4f)'%(mean(av[:3]), mean(av[:3])/Ntot))
plt.plot(ref_cutoff[:,0], ref_cutoff[:,1], 'r-', label=r'maximum criterion, $\delta t_c\approx 0.2$')
plt.xticks(delta_t)
plt.grid(axis='x')
plt.axis([0, 1.01, 1500, 2000])
plt.legend(loc = 'upper right', numpoints=1)
plt.xlabel('dimensionless topological time step (=topological time step / dissociation time)')
plt.ylabel('NAS')
plt.show()

# plt.clf()
# plt.ion()

# dt = 0.001
# Nw = 1000
# tb = dt*Nw

# ax = plt.subplot(111)

# for i in range(Nf):
#     ax.plot(dat[i][:,0]*tb/RT[i], dat[i][:,1], '-', color=csel[i], label='av=%4.2f(%6.4f), RT=%d'%(av[i], av[i]/float(Ntot), RT[i]))
# ax.legend(loc = 'lower right', prop=fp)
# ax.axis([0, 3, 0, 2000])
# ax.set_xlabel('topological time')
# ax.set_ylabel('NAS')
# ax2 = ax.twinx()
# ax2.grid()
# ax2.axis([0, 3, 0, 2000/Ntot])
# ax2.set_ylabel('fraction')

# plt.show()
