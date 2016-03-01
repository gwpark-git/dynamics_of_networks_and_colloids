


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




# fn_base = ('../NP0400_LD10P3_C100_SF15_RT100_longer/NP0400_LD10P3_C100')
# Np = 400
# N_dimension = 3
# box_dimension = 10.0
# lc = 0.12
# alpha = 1.5

if size(sys.argv) < 8:
    print 'USAGE:'
    print 'argv[1] == base file name'
    print 'argv[2] == out file name'
    print 'argv[3] == Np'
    print 'argv[4] == N_dimension'
    print 'argv[5] == box_dimension'
    print 'argv[6] == lc'
    print 'argv[7] == alpha'
    print 'argv[8] == minimum_time_index (DEFAULT 0)'
    print 'argv[9] == maximum_time_index (DEFAULT -1)'
else:

    fn_base = sys.argv[1]
    Np = int(sys.argv[3])
    N_dimension = int(sys.argv[4])
    box_dimension = float(sys.argv[5])
    lc = float(sys.argv[6])
    alpha = float(sys.argv[7])
    try:
        minimum_time_index = int(sys.argv[8])
        maximum_time_index = int(sys.argv[9])
    except:
        minimum_time_index = 0
        maximum_time_index = -1

    with open (fn_base + '.traj', 'r') as f_traj:
        with open (fn_base + '.hash', 'r') as f_index:
            with open (fn_base + '.weight', 'r') as f_weight:
                IDENTIFIER = 1
                # beta = []
                # r2_av = []
                dat = []
                t_cnt = 0
                while(IDENTIFIER):
                    # IDENTIFIER = 0
                    try:
                        if (t_cnt >= minimum_time_index and (t_cnt < maximum_time_index or maximum_time_index == -1)):
                            if (t_cnt % 100 == 0 and t_cnt <> 0):
                                print t_cnt, dat[-1][0], dat[-1][1]
                            t_cnt += 1
                            pos = read_traj_step(f_traj, Np, N_dimension)
                            connectivity = read_connectivity_step(f_index, f_weight, Np)
                            tmp_beta = 0.
                            tmp_cnt = 0
                            tmp_cnt_r2 = 0
                            tmp_r2 = 0.
                            for i in range(Np):
                                tmp_beta += connectivity[i, i]*1 # contribution for roof chain
                                tmp_cnt += connectivity[i, i]
                                for j in range(i + 1, Np):
                                    if connectivity[i, j] > 0:
                                        tmp_R = lin.norm(mapped_rel_vec_Rij(pos[i,:], pos[j,:], box_dimension))
                                        tmp_r2 += connectivity[i,j]*tmp_R**2.0
                                        tmp_cnt_r2 += connectivity[i,j]
                                        tmp_beta += connectivity[i,j]*compute_beta(tmp_R, N_dimension, lc, alpha)
                                        tmp_cnt += connectivity[i,j]
                            # beta.append(tmp_beta/float(tmp_cnt))
                            dat.append([tmp_r2/float(tmp_cnt_r2), tmp_beta/float(tmp_cnt)])
                        elif (t_cnt >= maximum_time_index):
                            break
                        else:
                            # this is set when the minimum accounting time index
                            pos = read_traj_step(f_traj, Np, N_dimension)
                            connectivity = read_connectivity_step(f_index, f_weight, Np)
                    except:
                        print 'END'
                        # beta = asarray(beta)
                        dat = asarray(dat)
                        break

    savetxt(sys.argv[2], dat)                            

