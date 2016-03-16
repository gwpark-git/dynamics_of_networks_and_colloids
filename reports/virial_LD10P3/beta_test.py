


from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as lin
import sys
sys.path.append('../post_processing')
from read_file import *


def F_Gaussian(r, N_dimension, alpha):
    return (N_dimension*alpha**2.0)*r

def compute_beta(r, N_dimension, lc, alpha):
    return exp(F_Gaussian(r, N_dimension, alpha)*lc)


fn_base = ('../NP0400_LD10P3_C100_SF15_RT100_longer/NP0400_LD10P3_C100')
Np = 400
N_dimension = 3
box_dimension = 10.0
lc = 0.12
alpha = 1.5
with open (fn_base + '.traj', 'r') as f_traj:
    with open (fn_base + '.hash', 'r') as f_index:
        with open (fn_base + '.weight', 'r') as f_weight:
            IDENTIFIER = 1
            beta = []
            t_cnt = 0
            while(t_cnt < 1000):
                # IDENTIFIER = 0
                if (t_cnt % 100 == 0):
                    print t_cnt, beta[-1]
                t_cnt += 1
                try:
                    pos = read_traj_step(f_traj, Np, N_dimension)
                    connectivity = read_connectivity_step(f_index, f_weight, Np)
                    tmp_beta = 0.
                    tmp_cnt = 0
                    for i in range(Np):
                        tmp_beta += connectivity[i, i]*1 # contribution for roof chain
                        tmp_cnt += connectivity[i, i]
                        for j in range(i + 1, Np):
                            if connectivity[i, j] > 0:
                                tmp_R = mapped_rel_vec_Rij(pos[i,:], pos[j,:], box_dimension)
                                tmp_beta += connectivity[i,j]*compute_beta(lin.norm(tmp_R), N_dimension, lc, alpha)
                                tmp_cnt += connectivity[i,j]
                    beta.append(tmp_beta/float(tmp_cnt))
                except:
                    print 'END'
                    beta = asarray(beta)
                    break
                
                            
                
