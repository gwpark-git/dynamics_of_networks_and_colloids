
from numpy import *
import matplotlib.pyplot as plt

# dt = dtype([('NC', float64), ('NAS', float64), ('f', float64), ])
           # 0: NC, 1:NAS, 2:f,   3:CN1, 4:CN2, 5:CN3, 6:RCN1, 7:RCN2, 8:RCN3, 9:R1, 10:R2, 11:R3
idx_NC   = 0
idx_NAS  = 1
idx_f    = 2
idx_CN1  = 3
idx_CN2  = 4
idx_CN3  = 5
idx_RCN1 = 6
idx_RCN2 = 7
idx_RCN3 = 8
idx_R1   = 9
idx_R2   = 10
idx_R3   = 11
dat = asarray([[5,   846,   2.64,  3.35, 5.68,  7.63,  0.98,   1.0,    0.99,   0.99, 2.04,  3.18],
               [10, 1998,   6.24,  4.11, 6.50,  8.19,  1.44,   1.27,   1.03,   0.98, 1.93,  2.84],
               [15, 3283,   10.26, 4.39, 6.47,  7.35,  1.76,   1.46,   1.06,   0.97, 1.83,  2.72],
               [20, 4589,   14.34, 4.46, 6.36,  6.59,  1.92,   1.53,   1.03,   0.96, 1.73,  2.60],
               [25, 5970,   18.66, 4.33, 5.65,  5.33,  1.95,   1.53,   0.90,   0.94, 1.67,  2.52]])


leg=[]
plt.clf()
plt.ion()
ax = plt.subplot(111)
leg.append(ax.plot(dat[:,idx_NC], dat[:,idx_NAS], 'ko-', markerfacecolor='white', markeredgecolor='black', label = 'Number of associations')[0])
ax.legend(loc = 'upper left', numpoints=1)
ax.set_xlabel('Number of chains per bead')
ax.set_ylabel('Number of associations (NAS)')
ax.axis([0, 30, 0, 6400])
ax.set_xticks(dat[:,idx_NC])
ax.set_yticks(dat[:,idx_NAS])
ax.set_yticklabels(dat[:,idx_NAS], rotation=90)
for i in range(shape(dat)[0]):
    ax.annotate('f=%s'%(dat[i,idx_f]), xy=(dat[i,idx_NC], dat[i,idx_NAS]), xytext=(dat[i,idx_NC]+0.5, dat[i,idx_NAS]), size='x-small')
# ax.annotate(dat[:,2], xy=(dat[:,0], dat[:,1]))
# ax.set_ticklabel(dat[:,1])
# ax.ticks_label(dat[:,1])
ax.grid(which='major', axis='x')

ax2 = ax.twinx()
ax2.plot(dat[:,idx_NC], dat[:, idx_RCN1], 'bs-', markerfacecolor='white', markeredgecolor='blue', label='RCN1')
for i in range(shape(dat)[0]):
    ax2.annotate('%4.2f (%4.2f)'%(dat[i, idx_RCN1], dat[i,idx_CN1]), xy=(dat[i,idx_NC], dat[i,idx_RCN1]), xytext=(dat[i,idx_NC]-0.01, dat[i,idx_RCN1]+0.01), color='blue', size='x-small')
ax2.plot(dat[:,0], dat[:, idx_RCN2], 'r^-', markerfacecolor='white', markeredgecolor='red', label='RCN2')
for i in range(shape(dat)[0]):
    ax2.annotate('%4.2f (%4.2f)'%(dat[i, idx_RCN2], dat[i,idx_CN2]), xy=(dat[i,idx_NC], dat[i,idx_RCN2]), xytext=(dat[i,idx_NC]-0.01, dat[i,idx_RCN2]+0.01), color='red', size='x-small')
ax2.plot(dat[:,0], dat[:, idx_RCN3], 'gd-', markerfacecolor='white', markeredgecolor='green', label='RCN3')
for i in range(shape(dat)[0]):
    ax2.annotate('%4.2f (%4.2f)'%(dat[i, idx_RCN3], dat[i,idx_CN3]), xy=(dat[i,idx_NC], dat[i,idx_RCN3]), xytext=(dat[i,idx_NC]-0.01, dat[i,idx_RCN3]+0.01), color='green', size='x-small')

ax2.legend(loc = 'lower right', numpoints=1)
ax2.axis([0, 30, 0.8, 2.0])
# ax2.axis([0, 30, 0, 6400*2/640])
# ax2.set_yticks(dat[:,2])
# ax2.set_yticklabels(dat[:,2], rotation=270)
# ax2.set_ylabel('functionality, 2*NAS/Np')

# plt.setp(ax, xticks=dat[:,0], xticklabels=dat[:,0])
# plt.setp(ax, yticks=dat[:,1], yticklabels=dat[:,1])
# plt.setp(ax2, yticks=dat[:,2], yticklabels=dat[:,2])
plt.show()

