from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../post_processing')
from read_file import *

def int_Gaussian(R_vec, alpha, N_dimension):
    return float(N_dimension)*(alpha**2.0)*R_vec

base_file_name = '../NP0400_LD10P3_C100'
Np = 400
box_dimension = 10.0
N_dimension = 3
alpha = 1.5
RF = []

initial_cnt = 2000

with open (base_file_name + '.traj', 'r') as f_traj:
    with open (base_file_name + '.hash', 'r') as f_index:
        with open (base_file_name + '.weight', 'r') as f_weight:
            IDENTIFIER = 1
            cnt = 0
            while(IDENTIFIER):
                try:
                    # IDENTIFIER = 0
                    cnt += 1
                    if (cnt > initial_cnt):
                        pos = read_traj_step(f_traj, Np, N_dimension)
                        connectivity = read_connectivity_step(f_index, f_weight, Np)
                        
                        # tmp_RF = zeros([6]) # for 3-dimensional case
                        RF_t = zeros([N_dimension, N_dimension])
                        for i in range(Np):
                            for j in range(i + 1, Np):
                                # remove inverse pairs
                                # note that it exclude j=i since force exerted on the loop chains are zero
                                # tmp_RF = zeros([N_dimension, N_dimension]) # temporal use
                                tmp_R = mapped_rel_vec_Rij(pos[i,:], pos[j,:], box_dimension)
                                tmp_F = int_Gaussian(tmp_R, alpha, N_dimension)
                                w_ij = connectivity[i, j]
                                for i_RF in range(N_dimension):
                                    for j_RF in range(N_dimension):
                                        RF_t[i_RF, j_RF] += w_ij * tmp_R[i_RF] * tmp_F[j_RF]
                        tmp_RF = RF_t.flatten()
                        RF.append(tmp_RF)
                        print cnt, tmp_RF
                    else:
                        tmp_pos = f_traj.readline()
                        for i in range(Np):
                            tmp_index = f_index.readline()
                            tmp_weight = f_weight.readline()
                except:
                    print 'end'
                    break
                

                            
                            
                                    

                            
                            
                    
    # and (open (base_file_name + '.hash', 'r') as f_index) and (open (base_file_name + '.weight', 'r') as f_weight)):

    
# pos = read_traj(base_file_name + '.traj', Np, N_dimension)
# connect = read_connectivity(base_file_name + '.hash', base_file_name + '.weight', Np)

