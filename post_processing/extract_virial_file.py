from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../post_processing')
from read_file import *

def int_Gaussian(R_vec, alpha, N_dimension):
    return float(N_dimension)*(alpha**2.0)*R_vec


# argv[1] == base file name
# argv[2] == output file name
# argv[3] == Np
# argv[4] == box_dimension
# argv[5] == N_dimension
# argv[6] == alpha
# argv[7] == initial cut (default is 1)

if size(sys.argv) < 7:
    print 'USAGE:'
    print 'argv[1] == base file name'
    print 'argv[2] == output file name'
    print 'argv[3] == Np'
    print 'argv[4] == box_dimension'
    print 'argv[5] == N_dimension'
    print 'argv[6] == alpha'
    print 'argv[7] == initial cut (default is 1)'
else:

    base_file_name = sys.argv[1]
    Np = int(sys.argv[3])
    box_dimension = float(sys.argv[4])
    N_dimension = int(sys.argv[5])
    alpha = float(sys.argv[6])
    RF = []
    try:
        initial_cnt = int(sys.argv[7])
    except:
        initial_cnt = 0
    
    with open (base_file_name + '.traj', 'r') as f_traj:
        with open (base_file_name + '.hash', 'r') as f_index:
            with open (base_file_name + '.weight', 'r') as f_weight:
                IDENTIFIER = 1
                cnt = 0
                while(IDENTIFIER):
                    try:
                        # IDENTIFIER = 0
                        cnt += 1
                        tmp_RF = zeros([N_dimension, N_dimension])
                        if (cnt >= initial_cnt):
                            pos = read_traj_step(f_traj, Np, N_dimension)
                            connectivity = read_connectivity_step(f_index, f_weight, Np)
                            # for i in range(N_dimension*N_dimension):
                            #     tmp_RF[i] = 0.
                            for i in range(Np):
                                for j in range(i + 1, Np):
                                    if connectivity[i,j] > 0:
                                        tmp_R = mapped_rel_vec_Rij(pos[i,:], pos[j,:], box_dimension)
                                        tmp_F = int_Gaussian(tmp_R, alpha, N_dimension)
                                        tmp_RF += connectivity[i,j] * outer(tmp_R, tmp_F)
                                    # remove inverse pairs
                                    # note that it exclude j=i since force exerted on the loop chains are zero
                                    # w_ij = connectivity[i, j]
                                    # tmp_RF = w_ij * outer(tmp_R, tmp_F)
                                    # cnt_ij = 0
                                        # for i_RF in range(N_dimension):
                                        #     for j_RF in range(N_dimension):
                                        #         tmp_RF[i_RF,j_RF] += connectivity[i,j] * tmp_R[i_RF] * tmp_F[j_RF]
                                    #         tmp_RF[cnt_ij] += w_ij * tmp_R[i_RF] * tmp_F[j_RF]
                                    #         cnt_ij += 1
                            RF.append(tmp_RF.flatten())
                            if cnt%100 == 1:
                                print cnt, RF[-1][:]
                        else:
                            tmp_pos = f_traj.readline()
                            for i in range(Np):
                                tmp_index = f_index.readline()
                                tmp_weight = f_weight.readline()

                    except:
                        print 'end count:', cnt
                        break
    savetxt(sys.argv[2], RF)



        # and (open (base_file_name + '.hash', 'r') as f_index) and (open (base_file_name + '.weight', 'r') as f_weight)):


    # pos = read_traj(base_file_name + '.traj', Np, N_dimension)
    # connect = read_connectivity(base_file_name + '.hash', base_file_name + '.weight', Np)

