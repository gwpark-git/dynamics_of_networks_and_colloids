


from numpy import *
import matplotlib.pyplot as plt

dat = []
dist_NP0400 = loadtxt('dist_bridge_NP0400_NC25_RT100.dat')[:,1]
beta_NP0400 = mean(dist_NP0400[-100:])
dat.append([400, beta_NP0400, 1./beta_NP0400, 0.162])

dist_NP0600 = loadtxt('dist_bridge_NP0600_NC25_RT100.dat')[:,1]
beta_NP0600 = mean(dist_NP0600[-100:])
dat.append([600, beta_NP0600, 1./beta_NP0600, 0.182])

# dist_NP0700 = loadtxt('dist_bridge_NP0700_NC25_RT100.dat')[:,1]
# beta_NP0700 = mean(dist_NP0700[-100:])
# dat.append([700, beta_NP0700, 1./beta_NP0700])

dist_NP0800 = loadtxt('dist_bridge_NP0800_NC25_RT100.dat')[:,1]
beta_NP0800 = mean(dist_NP0800[-100:])
dat.append([800, beta_NP0800, 1./beta_NP0800, 0.183])

dat = asarray(dat)
plt.close()
plt.ion()

colP = ['blue', 'green', 'red']

plt.figure(figsize=(6,8))

# plt.plot(1./dat[:,1], dat[:,3], 'b-')
for i in range(shape(dat)[0]):
    plt.plot(1./dat[i, 1],dat[i,3], 'o', markersize=10, markerfacecolor=colP[i], label = 'Np=%d'%(dat[i,0]))
plt.legend(loc = 'lower right', numpoints=1)
txt = []
for i in range(3):
    txt.append('%4.3f\n(%d)'%(1./dat[i,1], dat[i,0]))
# plt.xticks(dat[:,0])
plt.xticks(1./dat[:,1], txt)
# plt.yticks(1./dat[:,1], ['%4.3f'% val for val in 1./dat[:,1]])
plt.xticks(dat[:,2])
plt.yticks(dat[:,3], ['%4.3f'%val for val in dat[:,2]])
plt.grid()
plt.xlabel('average dissociation time')
plt.ylabel(r'measured correlation time')

# plt.figure(figsize=(6,8))

# plt.plot(dat[:,0], 1./dat[:,1], 'b-')
# for i in range(shape(dat)[0]):
#     plt.plot(dat[i, 0], 1./dat[i, 1], 'o', markersize=10, markerfacecolor=colP[i], label = 'Np=%d'%(dat[i,0]))
# plt.legend(loc = 'upper right', numpoints=1)
# # txt = []
# # for i in range(4):
# #     txt.append('%4.3f\n(%d)'%(1./dat[i,1], dat[i,0]))
# plt.xticks(dat[:,0])
# plt.yticks(1./dat[:,1], ['%4.3f'% val for val in 1./dat[:,1]])
# plt.grid()
# plt.xlabel('number of particles')
# plt.xlabel(r'number of particles')
# plt.ylabel(r'average dissociation time, $\bar{\tau} = 1/\bar{\beta}$')
# # plt.ylabel(r'average detachment frequency, $\bar{\beta}$')

# # plt.plot(dat[:,0], 1./dat[:,1], 'bo-')
# # plt.xticks(dat[:,0])
# # plt.yticks(1./dat[:,1])

# # plt.plot(dat[:,1], dat[:,2], 'bo-')

plt.show()
