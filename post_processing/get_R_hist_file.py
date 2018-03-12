
from matplotlib.font_manager import FontProperties
from numpy import *
# import matplotlib.pyplot as plt
import sys
from math import ceil

def index_array(x, dx):
    # this index function is on the purpose of specialized functionality with negative coordinate
    # for instance, if we have array such as | ... | -1 | 0 | 1 | .... |
    # the first incremental, say dx, must be dx/2 otherwise index 0 account the intensity of 2dx
    # from the second incremental (both of positive and negative), the incremental identity must be dx

    # the array of index is used such an way: |-2 | -1 | 0 | 1 | 2 | 
    # which is identical to the | 0 | 1 | 2 | -2 | -1 | (equivalently | 0 | 1 | 2 | 3 | 4 | )
    # since python handle index of array -1 as the last component of the array (if N is dimension of array, N-1 is the last component of array which is identical to the -1)

    # the number of total array must be accounted this aspect
    # for instance, if box_dimension is given by 10, then the maximum axis dimension of connector vector is half of box_dimension, 5.
    # let assume that dx = 2, then the number of positive index arrays can be int((5 - dx/2(=1))/2) = 2, which must be identical to the number of negative index.
    # in consequence, number of total arrays becomes 1(accounting zero) + 2 (#positive arr.) + 2 (#negative arr.) = 5, which must be odd number.

    # once all the coordination are specified, it can be easily used as "hist_x[index_array(x, dx)] += 1" weather positive or negative value of x.
    
    x = float(x)
    dx = float(dx)
    x0 = dx/2.
    return int(sign(x)*(abs(x) + x0)/dx)

def get_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
    return kd[argmin(abs(kd))]+x

def map_minimum_image_Rj_from_Ri(Ri, Rj, box_dimension):
    N_dimension = size(Rj)
    minimum_image_Rj = zeros(N_dimension)
    for i in range(N_dimension):
        minimum_image_Rj[i] = get_minimum_distance_k_from_x(Ri[i], Rj[i], box_dimension)
    return minimum_image_Rj

def rel_vec_Rij(Ri, Rj):
    return Rj - Ri

def RR_over_beads(pos, connectivity, box_dimension):
    Np, N_dimension = shape(pos)
    RR = zeros([N_dimension, N_dimension])
    for i in range(Np):
        tmp_index = connectivity[i,i:].nonzero()[0] + i
        #note that nonzero returns tuple. So, the [0] indice will return ndarray
        for j in tmp_index:
            R_j = map_minimum_image_Rj_from_Ri(pos[i,:], pos[j,:], box_dimension)
            R_ij = rel_vec_Rij(pos[i,:], R_j)
            RR += connectivity[i,j]*outer(R_ij, R_ij)
    RR /= float(Np)
    return RR

def hist_R_over_beads(pos, connectivity, box_dimension, hist_R, N_dimension, dr):
    Np, N_dimension = shape(pos)
    # RR = zeros([N_dimension, N_dimension])
    count = 0
    for i in range(Np):
        tmp_index = connectivity[i,i:].nonzero()[0] + i
        #note that nonzero returns tuple. So, the [0] indice will return ndarray
        for j in tmp_index:
            R_j = map_minimum_image_Rj_from_Ri(pos[i,:], pos[j,:], box_dimension)
            R_ij = rel_vec_Rij(pos[i,:], R_j)
            for k in range(N_dimension):
                hist_R[k, index_arr(R_ij[k], dr)] += 1
            count += 1
    return count


# Def R_over_beads(pos, connectivity, box_dimension):
#     Np, N_dimension = shape(pos)
#     R = zeros(


def read_traj(f_traj, Np, N_dimension):
    tmp_str = f_traj.readline().split('\t')[:-1]
    pos = zeros([Np, N_dimension])
    cnt = 1 # note that tmp_str[0] is time
    for i in range(Np):
        for k in range(N_dimension):
            pos[i, k] = float(tmp_str[cnt])
            cnt += 1
        # cnt += 2 # this will ignore velocity information
        cnt += N_dimension
    return pos

def read_connectivity(f_index, f_weight, Np):
    connectivity = zeros([Np, Np])
    str_index_table = []
    str_weight_table = []
    for i in range(Np):
        str_index_table.append(f_index.readline().split('\t')[:-1])
        str_weight_table.append(f_weight.readline().split('\t')[:-1])
    N_cols = shape(str_index_table)[1]
    for i in range(Np):
        index_i = int(str_index_table[i][0])
        for j in range(N_cols):
            index_j = int(str_index_table[i][j])
            if index_j == -1:
                break
            else:
                connectivity[index_i, index_j] = int(str_weight_table[i][j])
    return connectivity
    
# def Np_lines_parsing_index(seq_lines, Np, N_dimension):
#     connectivity = zeros([Np, Np])
    
# def line_parsing_hash(

if size(sys.argv) < 6:
    print 'USAGE:'
    print 'argv[1] == base filename'
    print 'argv[2] == out filename'
    print 'argv[3] == number of particles'
    print 'argv[4] == box dimension'
    print 'argv[5] == number of dimension (current version has symmetric dimensions)'
    print 'argv[6] == dx_k (current version has symmetric dimensions)'
    print 'argv[7] == Number of initial cuts (for getting steady-state data)'
    # print 'argv[8] == strider (to remove serially correlated data)'
else:
    fn_base = sys.argv[1]
    fn_traj = fn_base + '.traj'
    fn_index = fn_base + '.hash'
    fn_weight = fn_base + '.weight'
    fn_out = fn_base + '_hist_R.dat'
    Np = int(sys.argv[3])
    box_dimension = float(sys.argv[4])
    N_dimension = int(sys.argv[5])
    N_cols = 2*N_dimension*Np + 1
    tn = []

    N_cuts = int(sys.argv[7])
    # N_strider = int(sys.argv[8])

    # read directional increments
    dr = float(sys.argv[6])

    # get number of arrays for each directions
    # this is of important for the following procedure
    N_direction = int((box_dimension - dr/2.)/2.)
    Nr = 1 + 2*N_direction
    print 'initialized with %d-dimensional case with dr=%f (Nr=%d)'%(N_dimension, dr, Nr)
    # making connector vecotr histogram based on Cartesian coordinate
    # the first array index N_dimension refers the related direction, i.e., hist_R[0] is for x-axis since index 0 refers x
    # the second array index Nr is following the description of index_array function
    hist_R = zeros([N_dimension, Nr])
    

    
    # RR_hist = [] 
    # RR = zeros([N_dimension, N_dimension])
    with open(fn_traj, 'r') as f_traj:
        with open(fn_index, 'r') as f_index:
            with open(fn_weight, 'r') as f_weight:
                cnt = 0
                cnt_lines = 0
                while(cnt_lines >= N_cuts):
                    try:
                        pos = read_traj(f_traj, Np, N_dimension)
                        connectivity = read_connectivity(f_index, f_weight, Np)
                        # RR_t = RR_over_beads(pos, connectivity, box_dimension)
                        cnt += hist_R_over_beads(pos, connectivity, box_dimension, hist_R, N_dimension, dr)
                        # tn.append(cnt)
                        # RR_hist.append(RR_t.flatten()) # flatten transform 1d from Nd
                        # RR += RR_t
                        # cnt += 1
                        # for i in range(Np):
                        #     for j in range(Np):
                        #         tmp_RR = zeros([N_dimension, N_dimension])                                          for p in range(N_dimension):
                        #             for q in range(N_dimension):
                        #                 tmp_RR[p,q] = 
                    except:
                        # RR_hist = asarray(RR_hist)
                        # RR /= float(cnt)
                        # tn = asarray(tn)
                        break
    # dat = zeros([size(tn), N_dimension**2.0 + 1])
    # dat[:,0] = tn
    # dat[:, 1:] = RR
    dat = hist_R
    savetxt(fn_out, dat)


    # if N_dimension == 2:
    #     # ref_xx = asarray([[0, RR[0, 0]],
    #     #                   [tn[-1], RR[0, 0]]])
    #     # ref_xy = asarray([[0, RR[0, 1]],
    #     #                   [tn[-1], RR[0, 1]]])
    #     # ref_yy = asarray([[0, RR[1, 1]],
    #     #                   [tn[-1], RR[1, 1]]])
    #     # # plt.clf()
    #     # # plt.ion()
    #     # plt.plot(tn, RR_hist[:,0], 'b-', label = r'$\langle R_xR_x\rangle_{p}$')
    #     # plt.plot(tn, RR_hist[:,3], 'r-', label = r'$\langle R_yR_y\rangle_{p}$')
    #     # plt.plot(tn, RR_hist[:,1], 'g-', label = r'$\langle R_xR_y\rangle_{p}$')
    #     # plt.plot(ref_xx[:,0], ref_xx[:,1], 'b:', linewidth=3, label = r'$\langle R_xR_x\rangle_{p, t}$ = %6.2f'%(RR[0,0]))
    #     # plt.plot(ref_yy[:,0], ref_yy[:,1], 'r:', linewidth=3, label = r'$\langle R_yR_y\rangle_{p, t}$ = %6.2f'%(RR[1,1]))
    #     # plt.plot(ref_xy[:,0], ref_xy[:,1], 'g:', linewidth=3, label = r'$\langle R_xR_y\rangle_{p, t}$ = %6.2f'%(RR[0,1]))

    #     # font_prop = FontProperties()
    #     # font_prop.set_size('small')

    #     # plt.legend(loc = 'center left', ncol=2, prop=font_prop)
    #     # plt.xlabel('time steps')
    #     # plt.ylabel(r'$\mathbf{RR}$ averaged over all connection pairs, $\langle \mathbf{R}\mathbf{R}\rangle_{p}$')
    #     # plt.grid('on')
    #     # plt.savefig(fn_out)
    #     # plt.show()

    # else:
    #     ref_xx = asarray([[0, RR[0, 0]],
    #                       [tn[-1], RR[0, 0]]])
    #     ref_xy = asarray([[0, RR[0, 1]],
    #                       [tn[-1], RR[0, 1]]])
    #     ref_xz = asarray([[0, RR[0, 2]],
    #                       [tn[-1], RR[0, 2]]])
    #     ref_yy = asarray([[0, RR[1, 1]],
    #                       [tn[-1], RR[1, 1]]])
    #     ref_yz = asarray([[0, RR[1, 2]],
    #                       [tn[-1], RR[1, 2]]])
    #     ref_zz = asarray([[0, RR[2, 2]],
    #                       [tn[-1], RR[2, 2]]])
    #     plt.plot(tn, RR_hist[:,0], 'b-', label = r'$\langle R_xR_x\rangle_{p}$')
    #     plt.plot(tn, RR_hist[:,4], 'r-', label = r'$\langle R_yR_y\rangle_{p}$')
    #     plt.plot(tn, RR_hist[:,8], 'g-', label = r'$\langle R_zR_z\rangle_{p}$')
    #     plt.plot(tn, RR_hist[:,1], '-', color='cyan', label = r'$\langle R_xR_y \rangle_{p}$')
    #     plt.plot(tn, RR_hist[:,2], '-', color='purple', label = r'$\langle R_xR_z \rangle_{p}$')
    #     plt.plot(tn, RR_hist[:,5], '-', color='yellow', label = r'$\langle R_yR_z \rangle_{p}$')
    #     # plt.plot(tn, RR_hist[:,3], 'r-', label = r'$\langle R_yR_y\rangle_{p}$')
    #     # plt.plot(tn, RR_hist[:,1], 'g-', label = r'$\langle R_xR_y\rangle_{p}$')
    #     plt.plot(ref_xx[:,0], ref_xx[:,1], 'b:', linewidth=3, label = r'$\langle R_xR_x\rangle_{p, t}$ = %6.4f'%(RR[0,0]))
    #     plt.plot(ref_yy[:,0], ref_yy[:,1], 'r:', linewidth=3, label = r'$\langle R_yR_y\rangle_{p, t}$ = %6.4f'%(RR[1,1]))
    #     plt.plot(ref_zz[:,0], ref_zz[:,1], 'g:', linewidth=3, label = r'$\langle R_zR_z\rangle_{p, t}$ = %6.4f'%(RR[2,2]))
    #     plt.plot(ref_xx[:,0], ref_xy[:,1], ':', color='cyan', linewidth=3, label = r'$\langle R_xR_y\rangle_{p, t}$ = %6.4f'%(RR[0,1]))
    #     plt.plot(ref_yy[:,0], ref_xz[:,1], ':', color='purple', linewidth=3, label = r'$\langle R_xR_z\rangle_{p, t}$ = %6.4f'%(RR[0,2]))
    #     plt.plot(ref_zz[:,0], ref_yz[:,1], ':', color='yellow', linewidth=3, label = r'$\langle R_yR_z\rangle_{p, t}$ = %6.4f'%(RR[1,2]))

    #     font_prop = FontProperties()
    #     font_prop.set_size('small')

    #     plt.legend(loc = 'center right', ncol=2, prop=font_prop)
    #     plt.xlabel('time steps')
    #     plt.ylabel(r'$\mathbf{RR}$ averaged over all connection pairs, $\langle \mathbf{R}\mathbf{R}\rangle_{p}$')
    #     plt.grid('on')
    #     plt.savefig(fn_out)
        
