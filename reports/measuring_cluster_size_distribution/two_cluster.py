
from numpy import *
from cluster_dist import *

sc = []; so = []; IDPC=[]; IDPI=[]
Np = 40
Nd = 2
box_dimension = 10.0
hash = loadtxt('../../test_step/two_cluster.hash')
traj = loadtxt('../../test_step/two_cluster.traj')
pos = []
for k in range(Np):
    pos.append(traj[k*2*Nd + 1:k*2*Nd + 1 + Nd])
pos = asarray(pos)

size = DFS_percolation_edge_test(hash=hash, pos=pos, Ld=box_dimension, stack_component=sc, queue_order=so, index=0, IDPC=IDPC, IDPI=IDPI)
IDPC=asarray(IDPC)
IDPI=asarray(IDPI)
