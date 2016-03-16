
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
fp = FontProperties()
fp.set_size('x-small')

file_list = ['NP3200_LD20P3_C100_SF15_RT100.NAS',
             'NP1350_LD15P3_C100_SF15_RT100.NAS']

dat_3200 = loadtxt(file_list[0])
Ntot_3200 = 3200*25.
dat_1350 = loadtxt(file_list[1])
Ntot_1350 = 1350*25.

# RT = [100., 100.]
# Nf = size(file_list)
# dat = []
# for fn in file_list:
#     dat.append(loadtxt(fn))

# Ntot = 1350*25.
# av = zeros(Nf)
# for i in range(Nf):
#     dat[i][:, 0] -= dat[i][0, 0]
#     Ni = shape(dat[i])[0]
#     av[i] = mean(dat[i][Ni/2:,1])

csel = ['blue', 'cyan', 'red', 'purple', 'green', 'brown', 'black']



# plt.clf()
# ref_line = asarray([[0, mean(av[:3])],
#                     [1000, mean(av[:3])]])
                     
# plt.plot(RT, av, 'bo-', markerfacecolor='white', markeredgecolor='blue', label = 'NAS vs Rt')
# plt.plot(ref_line[:,0], ref_line[:,1], 'k-', label= 'ref line, %4.2f(%6.4f)'%(mean(av[:3]), mean(av[:3])/Ntot))
# plt.axis([0, 1000, 1500, 2000])
# plt.xticks(RT)
# plt.grid()
# plt.legend(loc = 'upper right', numpoints=1)
# plt.xlabel('Rational Number, topological time step/Brownian time step')
# plt.ylabel('NAS')
# plt.show()

plt.clf()
plt.ion()
dat_3200[:,0] -= dat_3200[0,0]
dat_1350[:,0] -= dat_1350[0,0]

# dt = 0.001
# Nw = 1000
# tb = dt*Nw
tb_3200 = 100*0.001
tb_1350 = 1000*0.001
# ax = plt.subplot(111)

# for i in range(Nf):
#     ax.plot(dat[i][:,0]*tb/RT[i], dat[i][:,1], '-', color=csel[i], label='av=%4.2f(%6.4f), RT=%d'%(av[i], av[i]/float(Ntot), RT[i]))
plt.plot(dat_3200[:,0]*tb/100., dat_3200[:,1]/Ntot_3200, 'b-', label='av=%4.2f(%6.4f), RT=%d, NP=%d'%(mean(dat_3200[shape(dat_3200)[0]/2:,1]), mean(dat_3200[shape(dat_3200)[0]/2:,1])/Ntot_3200, 100, 3200))
plt.plot(dat_1350[:,0]*tb/100., dat_1350[:,1]/Ntot_1350, 'r-', label='av=%4.2f(%6.4f), RT=%d, NP=%d'%(mean(dat_1350[shape(dat_1350)[0]/2:,1]), mean(dat_1350[shape(dat_1350)[0]/2:,1])/Ntot_1350, 100, 1350))
plt.legend(loc = 'lower right')
# plt.axis([0, 3, 0, 2000])
plt.xlabel('topological time')
plt.ylabel('fraction of NAS')
plt.grid()
plt.show()
