from numpy import *
import matplotlib.pyplot as plt

dat = []

# NAS0400 = loadtxt('../../NP0400_LD10P3_C100_SF15_RT100/NP0400_LD10P3_C100.NAS')[:,1]/(25*400.)
NP0400 = loadtxt('../../NP0400_LD10P3_C100_SF15_RT100/NP0400_LD10P3_C100.ener')
NAS0400 = NP0400[:,4]/(25*400.)
dat.append([400, mean(NAS0400[size(NAS0400)/2:])])

NP0532 = loadtxt('../../NP0532_LD11P3_C100_SF15_RT100/NP0532_LD11P3_C100.ener')
NAS0532 = NP0532[:,4]/(25*532.)
dat.append([532, mean(NAS0532[size(NAS0532)/2:])])

NP0691 = loadtxt('../../NP0691_LD12P3_C100_SF15_RT100/NP0691_LD12P3_C100.ener')
NAS0691 = NP0691[:,4]/(25*691.)
dat.append([691, mean(NAS0691[size(NAS0691)/2:])])


NP0879 = loadtxt('../../NP0879_LD13P3_C100_SF15_RT100/NP0879_LD13P3_C100.ener')
NAS0879 = NP0879[:,4]/(25*879.)
dat.append([879, mean(NAS0879[size(NAS0879)/2:])])

# NP1098 = loadtxt('../../NP1098_LD10P3_C100_SF15_RT100/NP1098_LD10P3_C100.ener')
# NAS1098 = NP1098[:,4]/532.

NAS1350 = loadtxt('../time_step_ratio/NP1350_LD15P3_C100_SF15_RT100.NAS')[:,1]/(25*1350.)
dat.append([1350, mean(NAS1350[size(NAS1350)/2:])])


dat = asarray(dat)
ref_line = asarray([[dat[0,0], mean(dat[1:,1])],
                    [dat[-1,0], mean(dat[1:,1])]])


plt.clf()
plt.ion()
plt.figure(figsize=(11,3))
plt.plot(dat[:,0], dat[:,1], 'bo-', markersize=10, markerfacecolor ='white', markeredgecolor='blue', label ='NAS faction')
plt.plot(ref_line[:,0], ref_line[:,1], 'k--', label='ref, %4.3f'%(mean(dat[1:,1])))
plt.xticks(dat[:,0], [10, 11, 12, 13, 15])
plt.legend(loc = 'lower right')
plt.grid(axis='x')
plt.axis([380, 1370, 0.044, 0.048])
plt.xlabel('dimensionless length for the box')
plt.ylabel('fraction of NAS')
plt.show()
# plt.plot(NAS0400, 'b-')
# plt.plot(NAS0532, 'r-')
# plt.plot(NAS0691, '
