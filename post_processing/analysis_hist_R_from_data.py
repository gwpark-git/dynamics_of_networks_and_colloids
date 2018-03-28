from numpy import *

from get_R_hist_file import *

def hist_R_over_beads_modified(pos, connectivity, box_dimension, hist_R, N_dimension, dr):
    # this is a modified version for the histogram analysis
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
                    hist_R[I_arr[0], I_arr[1], I_arr[2], 3] += connectivity[i, index_j_con]
                elif (N_dimension==2):
                    hist_R[I_arr[0], I_arr[1], 2] += connectivity[i, index_j_con]
                else:
                    print 'wrong dimensionality'

            count += 1
    return count



def measure_partition_fxyz(f_xyz, dx):
    # note that the given structure must be 4-dimensional array as 
    # shape(f_xyz) = [Nx, Ny, Nz, 4] where the last element (0,1,2) 
    # for coordinate and (3) for the real intensity
    
    # note that this function is not optimized currently
    # with optimization, we can reduce computational time by 2^3 where 3 is the dimensionality
    
    # note that currently only the symmetric increments is applied.
    Nx, Ny, Nz, NDp1 = shape(f_xyz)
    Z = 0.
    for i in range(Nx - 1):

        tmp_Zj_i0 = 0.
        tmp_Zj_i1 = 0.
        for j in range(Ny - 1):
            tmp_Zk_j0_i0 = 0.
            tmp_Zk_j1_i0 = 0.
        
            tmp_Zk_j0_i1 = 0.
            tmp_Zk_j1_i1 = 0.
            for k in range(Nz - 1):
                tmp_Zk_j0_i0 += 0.5*dx*(f_xyz[i, j, k+1, 3] + f_xyz[i,j,k, 3])
                tmp_Zk_j1_i0 += 0.5*dx*(f_xyz[i, j+1, k+1, 3] + f_xyz[i, j, k, 3])

                tmp_Zk_j0_i1 += 0.5*dx*(f_xyz[i+1, j, k+1, 3] + f_xyz[i+1,j,k, 3])
                tmp_Zk_j1_i1 += 0.5*dx*(f_xyz[i+1, j+1, k+1, 3] + f_xyz[i+1, j, k, 3])

            
            tmp_Zj_i0 += 0.5*dx*(tmp_Zk_j0_i0 + tmp_Zk_j1_i0)
            tmp_Zj_i1 += 0.5*dx*(tmp_Zk_j0_i1 + tmp_Zk_j1_i1)
        Z += 0.5*dx*(tmp_Zj_i0 + tmp_Zj_i1)
    return Z

def FENE_weight(q, q_max, alpha_factor):
    if (q < q_max):
        # q>= q_max is already prohibited by core of simulation code
        # the given array, however, is possibly to have q >= q_max, since it describes all the increment of spatial domain
        # which, however, the intensity always goes to zero.
        return 3.*alpha_factor**2.0 /(1. - (q/q_max)**2.0)
    return 0.

def FENE_weight_to_PDF(normalized_hist_R, q_max, alpha_factor):
    Nx, Ny, Nz, NDp1 = shape(normalized_hist_R)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                q = sqrt(normalized_hist_R[i,j,k, 0]**2.0 + normalized_hist_R[i,j,k,1]**2.0 + normalized_hist_R[i,j,k,2]**2.0)
                normalized_hist_R[i,j,k, 3] *= FENE_weight(q, q_max, alpha_factor)
    return 0

def measure_average_xy(normalized_f_XYZ):
    re = 0.
    f_XYZ = normalized_f_XYZ # just making new name
    Nx, Ny, Nz, NDp1 = shape(f_xyz)
    Z = 0.
    for i in range(Nx - 1):

        tmp_Zj_i0 = 0.
        tmp_Zj_i1 = 0.
        for j in range(Ny - 1):
            tmp_Zk_j0_i0 = 0.
            tmp_Zk_j1_i0 = 0.
        
            tmp_Zk_j0_i1 = 0.
            tmp_Zk_j1_i1 = 0.
            for k in range(Nz - 1):
                tmp_Zk_j0_i0 += 0.5*dx*(f_xyz[i, j, k+1, 3] + f_xyz[i,j,k, 3])
                tmp_Zk_j1_i0 += 0.5*dx*(f_xyz[i, j+1, k+1, 3] + f_xyz[i, j, k, 3])

                tmp_Zk_j0_i1 += 0.5*dx*(f_xyz[i+1, j, k+1, 3] + f_xyz[i+1,j,k, 3])
                tmp_Zk_j1_i1 += 0.5*dx*(f_xyz[i+1, j+1, k+1, 3] + f_xyz[i+1, j, k, 3])

            
            tmp_Zj_i0 += 0.5*dx*(tmp_Zk_j0_i0 + tmp_Zk_j1_i0)
            tmp_Zj_i1 += 0.5*dx*(tmp_Zk_j0_i1 + tmp_Zk_j1_i1)
        Z += 0.5*dx*(tmp_Zj_i0 + tmp_Zj_i1)
    return Z


def gen_hist_R_arr(Np, box_dimension, N_cuts, dx):
    N_dimension = 3
    N_cols = 2*N_dimension*Np + 1
    dr = dx
    N_direction = int((box_dimension/2. - dr/2.)/dr)
    Nr = 1 + 2*N_direction
    hist_R = zeros([Nr, Nr, Nr, N_dimension + 1])
    return hist_R

def get_intensity_R_from_data(fn_base, hist_R, Np, box_dimension, N_cuts, dx):
    fn_traj = fn_base + '.traj'
    fn_index = fn_base + '.hash'
    fn_weight = fn_base + '.weight'
    fn_out = fn_base + '_hist_R.dat'

    N_dimension = 3
    N_cols = 2*N_dimension*Np + 1
    tn = []


        # read directional increments
    dr = dx

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
#    hist_R = zeros([Nr, Nr, Nr, 4])

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
                            # if ((cnt_lines - N_cuts)%100 == 0):
                            #     print 'currently working with line number %d'%(cnt_lines)
                            cnt += hist_R_over_beads_modified(pos, connectivity, box_dimension, hist_R, N_dimension, dr)
                        cnt_lines += 1
                    except:
                        print '[break]ing line number = ', cnt_lines
                        break
        # note that the following codes are only compatible with 3-dimenional space

    return cnt_lines


def get_hist_R_from_list(fn_list, Np, box_dimension, N_cuts, dx):
    hist_R = gen_hist_R_arr(Np, box_dimension, N_cuts, dx)
    N_t = 0
    with open (fn_list, 'r') as f_list:
        print 'starting with ', fn_list
        count_samples = 0.
        time_stamp_long = 0
        tmp_time_stamp_long = 0

        for line in f_list:
            tmp_hist_R = copy(hist_R)            
            fn_base_tmp = line.replace('\n', '').replace('.ener','') # removing end-line deliminater and specific file name
            print '   currently working with ', fn_base_tmp.split('/')[-1]
            tmp_time_stamp_long = get_intensity_R_from_data(fn_base_tmp, hist_R, Np, box_dimension, N_cuts, dx)
            if (count_samples == 0):
                time_stamp_long = tmp_time_stamp_long
                Nt = time_stamp_long
                count_samples += 1.
                
            else:
                if tmp_time_stamp_long <> time_stamp_long: # when less or more time stamp of a condition is used
                    hist_R = copy(tmp_hist_R) 
                    print 'Caution: condition ', fn_base_tmp.split('/')[-1], 'has different length of time stamp. The histogram is excluded.'
                else: # correct one
                    count_samples += 1.
        
    hist_R[:, :, :, 3] /= float(count_samples)*float(Nt) # correcting intensity with number of samples and number of time
    return hist_R


if __name__ == "__main__":
    if size(sys.argv) < 6:
        print 'USAGE:'
        print 'argv[1] == filename of list'
        print 'argv[2] == out filename'
        print 'argv[3] == number of particles'
        print 'argv[4] == box dimension'
        print 'argv[5] == dx'
        print 'argv[6] == number of initial cuts'
    else:
        fn_list = sys.argv[1]
        fn_out = sys.argv[2]
        Np = int(sys.argv[3])
        box_dimension = float(sys.argv[4])
        dx = float(sys.argv[5])
        N_cuts = int(sys.argv[6])

        hist_R = get_hist_R_from_list(fn_list, Np, box_dimension, N_cuts, dx)
        save(fn_out, hist_R)
