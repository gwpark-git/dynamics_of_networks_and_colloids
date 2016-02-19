
from matplotlib.font_manager import FontProperties
from numpy import *
import matplotlib.pyplot as plt
import sys

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

if size(sys.argv) < 4:
    print 'USAGE:'
    print 'argv[1] == base filename'
    print 'argv[2] == number of particles'
    print 'argv[3] == box dimension'
    print 'argv[4] == number of dimension'
else:
    fn_base = sys.argv[1]
    fn_traj = fn_base + '.traj'
    fn_index = fn_base + '.hash'
    fn_weight = fn_base + '.weight'
    fn_out = fn_base + '_averaged_RR.pdf'
    Np = int(sys.argv[2])
    box_dimension = float(sys.argv[3])
    N_dimension = int(sys.argv[4])
    N_cols = 2*N_dimension*Np + 1
    tn = []
    RR_hist = [] 
    RR = zeros([N_dimension, N_dimension])
    with open(fn_traj, 'r') as f_traj:
        with open(fn_index, 'r') as f_index:
            with open(fn_weight, 'r') as f_weight:
                cnt = 0
                while(1):
                    try:
                        pos = read_traj(f_traj, Np, N_dimension)
                        connectivity = read_connectivity(f_index, f_weight, Np)
                        RR_t = RR_over_beads(pos, connectivity, box_dimension)
                        tn.append(cnt)
                        RR_hist.append(RR_t.flatten()) # flatten transform 1d from Nd
                        RR += RR_t
                        cnt += 1
                        # for i in range(Np):
                        #     for j in range(Np):
                        #         tmp_RR = zeros([N_dimension, N_dimension])                                          for p in range(N_dimension):
                        #             for q in range(N_dimension):
                        #                 tmp_RR[p,q] = 
                    except:
                        RR_hist = asarray(RR_hist)
                        RR /= float(cnt)
                        tn = asarray(tn)
                        break


    if N_dimension == 2:
        ref_xx = asarray([[0, RR[0, 0]],
                          [tn[-1], RR[0, 0]]])
        ref_xy = asarray([[0, RR[0, 1]],
                          [tn[-1], RR[0, 1]]])
        ref_yy = asarray([[0, RR[1, 1]],
                          [tn[-1], RR[1, 1]]])
        # plt.clf()
        # plt.ion()
        plt.plot(tn, RR_hist[:,0], 'b-', label = r'$\langle R_xR_x\rangle_{p}$')
        plt.plot(tn, RR_hist[:,3], 'r-', label = r'$\langle R_yR_y\rangle_{p}$')
        plt.plot(tn, RR_hist[:,1], 'g-', label = r'$\langle R_xR_y\rangle_{p}$')
        plt.plot(ref_xx[:,0], ref_xx[:,1], 'b:', linewidth=3, label = r'$\langle R_xR_x\rangle_{p, t}$ = %6.2f'%(RR[0,0]))
        plt.plot(ref_yy[:,0], ref_yy[:,1], 'r:', linewidth=3, label = r'$\langle R_yR_y\rangle_{p, t}$ = %6.2f'%(RR[1,1]))
        plt.plot(ref_xy[:,0], ref_xy[:,1], 'g:', linewidth=3, label = r'$\langle R_xR_y\rangle_{p, t}$ = %6.2f'%(RR[0,1]))

        font_prop = FontProperties()
        font_prop.set_size('small')

        plt.legend(loc = 'center left', ncol=2, prop=font_prop)
        plt.xlabel('time steps')
        plt.ylabel(r'$\mathbf{RR}$ averaged over all connection pairs, $\langle \mathbf{R}\mathbf{R}\rangle_{p}$')
        plt.grid('on')
        plt.savefig(fn_out)
        # plt.show()

    else:
        ref_xx = asarray([[0, RR[0, 0]],
                          [tn[-1], RR[0, 0]]])
        ref_xy = asarray([[0, RR[0, 1]],
                          [tn[-1], RR[0, 1]]])
        ref_xz = asarray([[0, RR[0, 2]],
                          [tn[-1], RR[0, 2]]])
        ref_yy = asarray([[0, RR[1, 1]],
                          [tn[-1], RR[1, 1]]])
        ref_yz = asarray([[0, RR[1, 2]],
                          [tn[-1], RR[1, 2]]])
        ref_zz = asarray([[0, RR[2, 2]],
                          [tn[-1], RR[2, 2]]])
        plt.plot(tn, RR_hist[:,0], 'b-', label = r'$\langle R_xR_x\rangle_{p}$')
        plt.plot(tn, RR_hist[:,4], 'r-', label = r'$\langle R_yR_y\rangle_{p}$')
        plt.plot(tn, RR_hist[:,8], 'g-', label = r'$\langle R_zR_z\rangle_{p}$')
        plt.plot(tn, RR_hist[:,1], '-', color='cyan', label = r'$\langle R_xR_y \rangle_{p}$')
        plt.plot(tn, RR_hist[:,2], '-', color='purple', label = r'$\langle R_xR_z \rangle_{p}$')
        plt.plot(tn, RR_hist[:,5], '-', color='yellow', label = r'$\langle R_yR_z \rangle_{p}$')
        # plt.plot(tn, RR_hist[:,3], 'r-', label = r'$\langle R_yR_y\rangle_{p}$')
        # plt.plot(tn, RR_hist[:,1], 'g-', label = r'$\langle R_xR_y\rangle_{p}$')
        plt.plot(ref_xx[:,0], ref_xx[:,1], 'b:', linewidth=3, label = r'$\langle R_xR_x\rangle_{p, t}$ = %6.4f'%(RR[0,0]))
        plt.plot(ref_yy[:,0], ref_yy[:,1], 'r:', linewidth=3, label = r'$\langle R_yR_y\rangle_{p, t}$ = %6.4f'%(RR[1,1]))
        plt.plot(ref_zz[:,0], ref_zz[:,1], 'g:', linewidth=3, label = r'$\langle R_zR_z\rangle_{p, t}$ = %6.4f'%(RR[2,2]))
        plt.plot(ref_xx[:,0], ref_xy[:,1], ':', color='cyan', linewidth=3, label = r'$\langle R_xR_y\rangle_{p, t}$ = %6.4f'%(RR[0,1]))
        plt.plot(ref_yy[:,0], ref_xz[:,1], ':', color='purple', linewidth=3, label = r'$\langle R_xR_z\rangle_{p, t}$ = %6.4f'%(RR[0,2]))
        plt.plot(ref_zz[:,0], ref_yz[:,1], ':', color='yellow', linewidth=3, label = r'$\langle R_yR_z\rangle_{p, t}$ = %6.4f'%(RR[1,2]))

        font_prop = FontProperties()
        font_prop.set_size('small')

        plt.legend(loc = 'center right', ncol=2, prop=font_prop)
        plt.xlabel('time steps')
        plt.ylabel(r'$\mathbf{RR}$ averaged over all connection pairs, $\langle \mathbf{R}\mathbf{R}\rangle_{p}$')
        plt.grid('on')
        plt.savefig(fn_out)
        
