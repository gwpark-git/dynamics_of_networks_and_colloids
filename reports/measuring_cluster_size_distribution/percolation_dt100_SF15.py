
from numpy import *
from cluster_dist import *
# sc = []; so = []; IDPC=[]; IDPI=[]

Np = 3200
Nd = 3
box_dimension = 20.0
hash = loadtxt('dt100.hash')
traj = loadtxt('dt100.traj')
pos = []
for k in range(Np):
    pos.append(traj[k*2*Nd + 1:k*2*Nd + 1 + Nd])
pos = asarray(pos)

# root_index = 1
# rc = []
# tmp_IDPC, tmp_IDPI = cluster_edge_DFS_travel_restricted_box_iter(hash = hash, pos = pos, Ld = box_dimension, record_component=rc, index = root_index)
# print_cluster_info(root_index, tmp_IDPC, tmp_IDPI, Nd)

# info = [0, 6, 10, 12, 18, 19, 22, 28, 40, 42, 44, 52, 55, 81, 100, 107, 115, 127, 144, 167, 175, 196, 223, 238, 243, 246, 268, 272, 282, 301, 307, 348, 353, 368, 379, 391, 415, 430, 433, 458, 461, 491, 539, 541, 565, 580, 632, 638, 761, 801, 979, 1129, 1159, 1225, 1435, 1875, 1929, 2305, 2693, 2941]

ident_dist = []
root_dist = []
size_dist = []
exist_queue = []
for i in range(Np):
# for i in info:
# for i in [1,2,3]:
    if i not in exist_queue:
        print 'root index = %d'%(i)
        ident_visit = 0
        rc = []
        IDPC, IDPI = cluster_edge_DFS_travel_restricted_box_iter(hash=hash, pos=pos, record_component = rc, Ld=box_dimension, index=i)
        
        for [j, k] in IDPI:
            if j not in exist_queue and k not in exist_queue:
                ident_visit = 1
                exist_queue.append(int(j))
                exist_queue.append(int(k))
        if ident_visit:
            root_dist.append(i)
            size_dist.append(measuring_cluster_size(IDPI))
            ident_dist.append(print_cluster_info(i, IDPC, IDPI, Nd, vervose=False))
            
        
Nsd = size(size_dist)
arg_map = argsort(size_dist)
dat = zeros([Nsd, 5])
per_x = []
per_y = []
per_z = []
# 0==root index, 1==size_dist, 2==per_x, 3==per_y, 4==per_z

for i in range(Nsd):
    dat[i, 0] = root_dist[arg_map[i]]
    dat[i, 1] = size_dist[arg_map[i]]
#     for k in range(Nd):
#         if ident_dist[arg_map[i]][0] > 0:
#             per_x.append([i, root_dist[arg_map[i]], size_dist[arg_map[i]]])
#         if ident_dist[arg_map[i]][1] > 0:
#             per_y.append([i, root_dist[arg_map[i]], size_dist[arg_map[i]]])
#         if ident_dist[arg_map[i]][2] > 0:
#             per_z.append([i, root_dist[arg_map[i]], size_dist[arg_map[i]]])

# per_x = asarray(per_x); per_y = asarray(per_y); per_z = asarray(per_z)                


savetxt('dt100_CD.dat', dat)    

    
# import matplotlib.pyplot as plt
# plt.clf()
# plt.ion()
# # plt.plot(root_dist, size_dist, 'b.-', label='cluster size dist.')
# plt.plot(dat[:,1], 'k.-', label = 'cluster size dist.')
# # plt.plot(per_x[:,0], per_x[:,2], 'ro', markersize=10, markerfacecolor='none', markeredgecolor='red', label = 'percolation along x')
# # plt.plot(per_y[:,0], per_y[:,2], 'b^', markersize=10, markerfacecolor='none', markeredgecolor='blue', label = 'percolation along y')
# # plt.plot(per_z[:,0], per_z[:,2], 'gv', markersize=10, markerfacecolor='none', markeredgecolor='green', label = 'percolation along z')

# plt.xlabel('index for distingushable cluster')
# plt.ylabel('size of cluster (= number of particles in the cluster)')
# plt.grid()
# plt.legend(loc = 'upper left', numpoints=1)
# plt.show()
