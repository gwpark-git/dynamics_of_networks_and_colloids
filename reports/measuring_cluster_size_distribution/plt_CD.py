
from numpy import *
import matplotlib.pyplot as plt

dat_dt050 = loadtxt('dt050_CD.dat')
N050 = size(dat_dt050[:,1])
dat_dt100 = loadtxt('dt100_CD.dat')
N100 = size(dat_dt100[:,1])

NC_max_dt050 = int(max(dat_dt050[:,1]))
NC_max_dt100 = int(max(dat_dt100[:,1]))
NC_max = max(NC_max_dt050, NC_max_dt100)
dt = NC_max/30.

# dat = zeros([N_max, 3])
# cnt_050 = 0; cnt_100 = 0;
P_dt050 = zeros(NC_max_dt050)
# for i in range(NC_max_dt050):
for i in range(size(dat_dt050[:,1])):
    index = int(dat_dt050[i, 1]/dt)
    P_dt050[index] += 1

P_dt100 = zeros(NC_max_dt100)
for i in range(size(dat_dt100[:,1])):
    index = int(dat_dt100[i, 1]/dt)
    P_dt100[index] += 1

    
    
plt.clf()
plt.ion()

ax = plt.subplot(111)
ax.plot(P_dt050, 'ro-', linewidth=1, markerfacecolor='white', markeredgecolor='red', label = 'Rt=50,  total = %d'%(sum(dat_dt050[:,1])))
ax.plot(P_dt100, 'bo-', linewidth=1, markerfacecolor='white', markeredgecolor='blue', label = 'Rt=100, total = %d'%(sum(dat_dt100[:,1])))

ax.grid()
ax.legend(loc = 'center right', numpoints=1)
ax.set_xlabel('cluster size')
ax.set_ylabel('number of clusters')
ax.axis([0, NC_max, -1, 350])
ax2 = ax.twinx()
ax2.plot(dat_dt050[:,1], range(size(dat_dt050[:,1])), 'r--', markersize=2, markerfacecolor='white', markeredgecolor='red')
ax2.plot(dat_dt100[:,1], range(size(dat_dt100[:,1])), 'b--', markersize=2, markerfacecolor='white', markeredgecolor='blue')
ax2.set_ylabel('sorted index for cluster')
ax2.axis([0, NC_max, 0, 900])


plt.show()
