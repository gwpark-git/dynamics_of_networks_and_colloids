import os, sys
lib_path = os.path.abspath(os.path.join('..','post_processing'))
sys.path.append(lib_path)
from lib_rdf import *
from scipy.linalg import norm
from numpy import *
import matplotlib.pyplot as plt

def pos(traj, t, i, k, Nd):
    return traj[t, 2*Nd*i + 1 + k]


if size(sys.argv) < 9:
    print 'USAGE:'
    print 'argv[1] == base file name'
    print 'argv[2] == output ddf file name'
    print 'argv[3] == starting time step'
    print 'argv[4] == stride'
    print 'argv[5] == end time step. default: -1 for end of data'
    print 'argv[6] == N_dimension'
    print 'argv[7] == number of particles'
    print 'argv[8] == box_dimension'
    print ' the following for asymmetric PBC box'
    print 'argv[9] == Lx'
    print 'argv[10] == Ly'
    print 'argv[11] == Lz'
else:
    fn_traj = sys.argv[1] + '.traj'
    fn_hash = sys.argv[1] + '.hash'
    fn_weight = sys.argv[1] + '.weight'
    fn_ddf_out = sys.argv[2]
    t_st = int(sys.argv[3])
    # N_t_block = int(sys.argv[4])
    N_stride = int(sys.argv[4])
    t_end = int(sys.argv[5])
    # fac_t = int(sys.argv[5])
    # dt = float(sys.argv[6])
    traj = loadtxt(fn_traj)
    if t_end == -1:
        t_end = shape(traj)[0] 
    # print t_st, t_end
    N_dimension = int(sys.argv[6])


    Np = int(sys.argv[7])

    box_dimension = zeros(N_dimension)
    for i in range(N_dimension):
        if size(sys.argv)<10:
            box_dimension[i] = float(sys.argv[8])
        else:
            box_dimension[i] = float(sys.argv[9+i])

    ddf = []

    # rho = 0.4
    cut_ratio = 0.5
    # fac_t = 1


    traj = loadtxt(fn_traj)
    N_cols = 0
    with open (fn_traj, 'r') as f_traj:
        with open (fn_hash, 'r') as f_hash:
            with open (fn_weight, 'r') as f_weight:
                for t in range(t_st, t_end, N_stride): #check definition for ts
                    # hash_st = t*Np # direction for initial position of hash, which is concide with weight.
                    tmp_ddf = [traj[t, 0]]
                    for i in range(Np):
                        # print t, t_end, i, tmp_hash
                        hash = map(int, f_hash.readline().replace('\t\n', '').split('\t'))

                        weight = map(int, f_weight.readline().replace('\t\n', '').split('\t'))
                        for hj, j in enumerate(hash):
                            if j>i: # prevent duplication bridges and itself (hj=0 => j=i)
                                rb_vec = get_rel_vec(traj, t, i, j, N_dimension, box_dimension) # asymmetric_safety
                                # ddf.append(float(weight[hj])*norm(rb_vec))
                                tmp_ddf.append(norm(rb_vec))
                    N_cols = max(size(tmp_ddf), N_cols)
                    ddf.append(tmp_ddf)
    re = zeros([shape(ddf)[0], N_cols])
    for i in range(shape(ddf)[0]):
        for j in range(size(ddf[i])):
            re[i,j] = ddf[i][j]
    savetxt(fn_ddf_out, re)
