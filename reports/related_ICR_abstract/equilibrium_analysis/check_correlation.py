


from numpy import *
import matplotlib.pyplot as plt

dat = asarray([[400, 1.026, 0.162],
               [600, 1.041, 0.182],
               [700, 1.046, 0.263],
               [800, 1.056, 0.183]])

plt.clf()
plt.ion()

colP = ['blue', 'green', 'brown', 'red']

plt.figure(figsize=(6,8))
plt.plot(1./dat[:,1], dat[:,2], 'b-')
for i in range(shape(dat)[0]):
    plt.plot(1./dat[i, 1], dat[i, 2], 'o', markersize=10, markerfacecolor=colP[i], label = 'Np=%d'%(dat[i,0]))
plt.legend(loc = 'upper right', numpoints=1)
txt = []
for i in range(4):
    txt.append('%4.3f\n(%d)'%(1./dat[i,1], dat[i,0]))
plt.xticks(1./dat[:,1], txt)
plt.yticks(dat[:,2], ['%4.3f'% val for val in dat[:,2]])
plt.grid()
# plt.xlabel('number of particles')
plt.xlabel(r'$\bar{\tau} = 1/\bar{\beta}$ (number of particles)')
plt.ylabel(r'correlation time')

# plt.plot(dat[:,0], 1./dat[:,1], 'bo-')
# plt.xticks(dat[:,0])
# plt.yticks(1./dat[:,1])

# plt.plot(dat[:,1], dat[:,2], 'bo-')

plt.show()
