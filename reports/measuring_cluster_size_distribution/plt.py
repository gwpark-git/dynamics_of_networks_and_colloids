
from numpy import *
import matplotlib.pyplot as plt

skey = argsort(size_dist)

Ndist = size(size_dist)
dat = zeros([Ndist, 2])
for i in range(size_dist):
    dat[i, 0] = root_dist[skey[i]]
    dat[i, 1] = size_dist[skey[i]]

plt.clf()
plt.ion()
plt.plot(dat[:, 1], 'b.', label = 'size distribution')
plt.ylabel('size of cluster')
plt.show()
    
