
from numpy import *
from cluster_dist import *
# sc = []; so = []; IDPC=[]; IDPI=[]

Np = 3200
Nd = 3
box_dimension = 20.0
hash = loadtxt('test_SF15.hash')
traj = loadtxt('test_SF15.traj')
pos = []
for k in range(Np):
    pos.append(traj[k*2*Nd + 1:k*2*Nd + 1 + Nd])
pos = asarray(pos)

# size = cd.cluster_edge_DFS_travel_restricted_box_iter(hash = hash, pos = pos, Ld = box_dimension, stack_component=sc, queue_order = so, index = 1, IDPC = IDPC, IDPI = IDPI)
# cd.print_cluster_info(root_index, IDPC, IDPI, Nd)

root_dist = []
size_dist = []
exist_queue = []
for i in range(10):
    if i not in exist_queue:
        ident_visit = 0
        sc = []; so = []; IDPC=[]; IDPI=[]
        DF_percolation_edge_restricted_box_iter(hash=hash, pos=pos, Ld=box_dimension, stack_component=sc, queue_order=so, index=i, IDPC=IDPC, IDPI=IDPI)
        cd.print_cluster_info(i, IDPC, IDPI, Nd)
        for [j, k] in IDPI:
            if j not in exist_queue:
                ident_visit = 1
                exist_queue.append(int(j))
            if k not in exist_queue:
                ident_visit = 1
                exist_queue.append(int(k))
        if ident_visit:
            root_dist.append(i)
            size_dist.append(shape(IDPC)[0])
        

