
from matplotlib.font_manager import FontProperties
from numpy import *
# import matplotlib.pyplot as plt
import sys
from math import ceil

def map_coordinate_from_index(index_arr, dr, N_dimension, N_direction):
    # Note that Nr represents number of increments for positive direction
    # total number of index are Nr = 1 + 2*N_direction

    # the remapping procedure strictly depends on the index handling that is described in index_array function.
    # here, the remapping procedure is way of inversion where the coordinate remap through an average between two incremental limits.
    # for instance, the first index (0) will remap to coordinate zero because (-dx/2 + dx/2)/2 = 0
    # the second index (1) will remap to dx since (dx/2 + (dx/2 + dx))/2 = dx.

    # in the case of negative index, we can identify the negative index through int(index/(N_direction + 1))==1 since from 1 to N_direction, we have N_direction-number of positive index,
    # then from N_direction + 1 to 2 N_direction, we have N_direction-number of negative index.
    # note that index for array is always positive while its original mapping (designed from negative to posistive index of array) should be positive integer
    
    dr = float(dr)
    r0 = dr/2.
    r_arr = zeros(N_dimension)
    for k in range(N_dimension):
        if int(index_arr[k]/(N_direction + 1)) == 1:
            tmp_index = index_arr[k] - N_direction
            r_arr[k] = -dr*(N_direction + 1 - tmp_index)
            # r_arr[k] = -dr*(index_arr[k]%(N_direction + 1))
        else:
            r_arr[k] = dr*index_arr[k]
    return r_arr
    

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

# def index_map_to_1d(I_arr, N_arr):
#     # I_arr = [Ix, Iy, Iz] (or [Ix, Iy] for 2d)
#     # N_arr = [Nx, Ny, Nz] which is number of arrays for each directions (or [Nx, Ny] for 2d)

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
    I_arr = [0, 0, 0] # integer array
    for i in range(Np):
        tmp_index = connectivity[i,i:].nonzero()[0] + i # this will have the index array
        #note that nonzero returns tuple. So, the [0] indice will return ndarray
        for cnt_index, j in enumerate(tmp_index): # (cnt_index + i + 1) is the index for column in connectivity matrix
            if j <> i:
                index_j_con = cnt_index + i + 1
                # excluding the roof chains
                # this explicit excluding will boost the performance of code dramatically since most of chains are in the roof status
                # in addition, the excluding is of important to distingush the intensity of (0, 0, 0) coordinate is only accounted for the bridge chains
                # which in generally almost zero because of excluded volume effect
                R_j = map_minimum_image_Rj_from_Ri(pos[i,:], pos[j,:], box_dimension)
                R_ij = rel_vec_Rij(pos[i,:], R_j)
                for k in range(N_dimension):
                    I_arr[k] = index_array(R_ij[k], dr)
                if (N_dimension==3):
                    hist_R[I_arr[0], I_arr[1], I_arr[2]] += connectivity[i, index_j_con]
                elif (N_dimension==2):
                    hist_R[I_arr[0], I_arr[1]] += connectivity[i, index_j_con]
                else:
                    print 'wrong dimensionality'
            
            count += 1
    return count


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
    # str_index_table = []
    # str_weight_table = []
    for i in range(Np):
        tmp_arr_index = f_index.readline().split('\t')[:-1]
        tmp_arr_weight = f_weight.readline().split('\t')[:-1]
        # str_index_table.append(tmp_arr_index)
        # str_weight_table.append(tmp_arr_weight)
        N_cols = size(tmp_arr_index)
        index_i = int(tmp_arr_index[0])
        for j in range(N_cols):
            index_j = int(tmp_arr_index[j])
            if index_j == -1:
                break
            else:
                connectivity[index_i, index_j] = int(tmp_arr_weight[j])
    return connectivity
    
if __name__ == "__main__":
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
        # note that box_dimension/2 is the maximum value of directional coordinate
        # because of cut-off scheme of the core simulation
        N_direction = int((box_dimension/2. - dr/2.)/dr)
        Nr = 1 + 2*N_direction
        print 'initialized with %d-dimensional case with dr=%f (Nr=%d)'%(N_dimension, dr, Nr)
        # index 0: middle x
        # index 1: middle y
        # index 2: middle z
        # index 3: intensity

        # making connector vecotr histogram based on Cartesian coordinate
        # the first array index N_dimension refers the related direction, i.e., hist_R[0] is for x-axis since index 0 refers x
        # the second array index Nr is following the description of index_array function
        # hist_R = zeros([N_dimension, Nr])
        hist_R = zeros([Nr, Nr, Nr])

        # RR_hist = [] 
        # RR = zeros([N_dimension, N_dimension])
        with open(fn_traj, 'r') as f_traj:
            with open(fn_index, 'r') as f_index:
                with open(fn_weight, 'r') as f_weight:
                    cnt = 0
                    cnt_lines = 0
                    if (cnt_lines < N_cuts):
                        for tmp_iter in xrange(N_cuts):
                            f_traj.readline()
                            cnt_lines += 1
                            for tmp_iter_2 in xrange(Np):
                                f_index.readline()
                                f_weight.readline()

                    while(1):
                        try:
                            pos = read_traj(f_traj, Np, N_dimension)
                            connectivity = read_connectivity(f_index, f_weight, Np)
                            # RR_t = RR_over_beads(pos, connectivity, box_dimension)
                            if(cnt_lines == N_cuts):
                                print 'line number %d meet the starting condition'%(N_cuts)
                            if(cnt_lines >= N_cuts):
                                if ((cnt_lines - N_cuts)%100 == 0):
                                    print 'currently working with line number %d'%(cnt_lines)
                                cnt += hist_R_over_beads(pos, connectivity, box_dimension, hist_R, N_dimension, dr)
                            cnt_lines += 1
                        except:
                            print 'breaking line number = ', cnt_lines
                            break
        # note that the following codes are only compatible with 3-dimenional space
        dat = []
        # print nonzero(hist_R)
        tmp_arr = zeros(4)
        cnt_line = 0
        for i in range(Nr):
            for j in range(Nr):
                for k in range(Nr):
                    if int(hist_R[i,j,k]) > 0:
                        # this if phrase to reduce the extracted data
                        # because hist_R will be large-sparse matrix
                        # most of the data can be ignored since the intensity is zero
                        tmp_arr[:3] = map_coordinate_from_index([i, j, k], dr, N_dimension, N_direction)
                        tmp_arr[3] = hist_R[i, j, k]
                        dat.append(tmp_arr.tolist())
                        # from numpy array to list function is important to prevent the change of data inside arrays
                        # note that dat was initialized with an empty list in order to use the member function .append
                        # this procedure has a benefit to compress the sparse matrix hist_R into dat
                        # then, the dat array will explicitly type change when we call savetxt function
                        cnt_line += 1
        # print dat
        print shape(dat)
        savetxt(fn_out, asarray(dat))


