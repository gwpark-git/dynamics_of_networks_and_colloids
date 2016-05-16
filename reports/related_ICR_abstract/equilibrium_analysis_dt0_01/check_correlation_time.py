


from numpy import *
import matplotlib.pyplot as plt

dat = []
dist_NP0400 = loadtxt('dist_NP0400_NC25_RT100.dat')[:,1]
beta_NP0400 = mean(dist_NP0400[-100:])
dat.append([400, beta_NP0400, 1./beta_NP0400, 0.221, 0.137])

dist_NP0600 = loadtxt('dist_NP0600_NC25_RT100.dat')[:,1]
beta_NP0600 = mean(dist_NP0600[-100:])
dat.append([600, beta_NP0600, 1./beta_NP0600, 0.190, 0.407])


dist_NP0800 = loadtxt('dist_NP0800_NC25_RT100.dat')[:,1]
beta_NP0800 = mean(dist_NP0800[-100:])
dat.append([800, beta_NP0800, 1./beta_NP0800, 0.214, 2.133])

dist_NP1000 = loadtxt('dist_NP1000_RT100_new.dat')[:,1]
beta_NP1000 = mean(dist_NP1000[-100:])
dat.append([1000, beta_NP1000, 1./beta_NP1000, 0.409, 2.687])

dist_NP1200 = loadtxt('dist_NP1200_RT100_new.dat')[:,1]
beta_NP1200 = mean(dist_NP1200[-100:])
dat.append([1200, beta_NP1200, 1./beta_NP1200, 0.422, 0.207])

dist_NP1400 = loadtxt('dist_NP1400_RT100_new.dat')[:,1]
beta_NP1400 = mean(dist_NP1400[-100:])
dat.append([1400, beta_NP1400, 1./beta_NP1400, 0.526, 1.114])


dat = asarray(dat)
plt.close()
plt.ion()

colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple']

plt.figure(figsize=(6,8))

# plt.plot(1./dat[:,1], dat[:,3], 'b-')
# for i in range(shape(dat)[0]):
#     plt.plot(1./dat[i, 1],dat[i,3], 'o', markersize=10, markerfacecolor=colP[i], label = 'Np=%d'%(dat[i,0]))
# plt.legend(loc = 'lower right', numpoints=1)
# txt = []
# for i in range(3):
#     txt.append('%4.3f\n(%d)'%(1./dat[i,1], dat[i,0]))
# # plt.xticks(dat[:,0])
# plt.xticks(1./dat[:,1], txt)
# # plt.yticks(1./dat[:,1], ['%4.3f'% val for val in 1./dat[:,1]])
# plt.xticks(dat[:,2])
# plt.yticks(dat[:,3], ['%4.3f'%val for val in dat[:,2]])
# plt.grid()
# plt.xlabel('average dissociation time')
# plt.ylabel(r'measured correlation time')


plt.plot(dat[:,0], dat[:,4], 'b-')
for i in range(shape(dat)[0]):
    plt.plot(dat[i,0], dat[i,4], 'o', markersize=10, markerfacecolor=colP[i], label = 'Np=%d, tc_3rd=%4.3f'%(dat[i,0], dat[i,4]))
plt.legend(loc = 'upper left', numpoints=1)
# txt = []
# for i in range(4):
#     txt.append('%4.3f\n(%d)'%(1./dat[i,1], dat[i,0]))
# plt.xticks(dat[:,0], ['%4.3f'% val for val in dat[:,0]])
plt.xticks(dat[:,0])
plt.yticks(dat[:,4], ['%4.3f'% val for val in dat[:,4]])
plt.grid()
plt.xlabel('number of particles')
# plt.xlabel(r'average dissociation time, $\bar{\tau}=1/\bar{\beta}$')
# plt.ylabel(r'average dissociation time, $\bar{\tau} = 1/\bar{\beta}$')
# plt.ylabel(r'average detachment frequency, $\bar{\beta}$')
plt.ylabel('correlation time')
plt.axis([380, 1420, 0.11, 4])
# plt.axis([0.448, 0.457, 0.18, 0.3])
# plt.axis([0.448, 0.457, 0.11, 3.0])
# plt.axis([380, 1020, 0, 3])

# plt.axis([380, 1020, 0.18, 0.3])
# # plt.plot(dat[:,0], 1./dat[:,1], 'bo-')
# # plt.xticks(dat[:,0])
# # plt.yticks(1./dat[:,1])

# # plt.plot(dat[:,1], dat[:,0], 'bo-')

plt.show()
