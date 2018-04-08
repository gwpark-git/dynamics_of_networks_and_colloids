from numpy import *
from analysis_hist_R_from_data import *
from scipy.linalg import norm


def coord_XYZ_to_SPE(XYZ_arr):
    r = norm(XYZ_arr)
    theta = arccos(XYZ_arr[2]/r)
    phi = arctan2(XYZ_arr[1], XYZ_arr[0])
    return asarray([r, theta, phi])

def index_array_SPE(val, d_val):
    # as described in other modified version in this script file, the data structure is differ from Cartesian one.
    # in consequence, the index functions has different definition where the half division is not necessary in this aspect
    # note again that we are keeping the same incremental values (say d_val) for all given val.
    # in addition, the given values are already sorted with respect to given index functions
    val = float(val)
    d_val = float(d_val)
    return int(val/d_val)


def hist_SPE_over_beads_modified(pos, connectivity, box_dimension, hist_SPE, N_div, M0):
    Np, N_dimension = shape(pos)
    count = 0
    I_arr = [0, 0, 0]
    dval_arr = asarray([hist_SPE[1,0,0,0] - hist_SPE[0,0,0,0], hist_SPE[0,1,0,1] - hist_SPE[0,0,0,1], hist_SPE[0,0,1,2] - hist_SPE[0,0,0,2]])

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
                # R_j = map_minimum_image_Rj_from_Ri(pos[i,:] , pos[j,:], box_dimension)
                R_j = map_minimum_image_Rj_from_Ri_simple_shear_3d(pos[i,:], pos[j,:], box_dimension, M0)
                R_ij = rel_vec_Rij(pos[i,:], R_j)
                R_ij_SPE = coord_XYZ_to_SPE(R_ij)
                tmp_r = norm(R_ij)
                if (tmp_r >= 3.0):
                    print 'Warning: given r is ', tmp_r
                    print '\t R_i =', pos[i,:], '\t Rj = ', pos[j, :]
                    print '\t R_ij = ', R_ij, ' new R_j = ', R_j, ', M0 = ', M0
                for k in range(N_dimension):
                    I_arr[k] = index_array_SPE(R_ij_SPE[k], dval_arr[k])
                if (N_dimension==3):
                    # hist_R[I_arr[0], I_arr[1], I_arr[2], 3] += connectivity[i, index_j_con]
                    # count += connectivity[i, index_j_con]
                    hist_R[I_arr[0], I_arr[1], I_arr[2], 3] += connectivity[i, j]
                    count += connectivity[i, j]
                    
                elif (N_dimension==2):
                    # hist_R[I_arr[0], I_arr[1], 2] += connectivity[i, index_j_con]
                    hist_R[I_arr[0], I_arr[1], 2] += connectivity[i, j]
                    
                else:
                    print 'wrong dimensionality'

    return count


def get_intensity_SPE_from_data(fn_base, hist_SPE, Np, box_dimension, N_cuts, N_div, Wi_R, Delta_t_strider):
    # the treatment of index and coordinate functions of spherical coordinate system has totally different nature
    # compared with Cartesian one. The main reason will be described in an individual description files for future
    # apart from the designed architecture, here, we will use a simplified version as a starting point.
    
    fn_traj = fn_base + '.traj'
    fn_index = fn_base + '.hash'
    fn_weight = fn_base + '.weight'
    fn_out = fn_base + '_hist_R.dat'

    N_dimension = 3
    tn = []
    
    with open(fn_traj, 'r') as f_traj:
        with open(fn_index, 'r') as f_index:
            with open(fn_weight, 'r') as f_weight:
                cnt_lines = 0
                if (cnt_lines < N_cuts):
                    # if cnt_lines > 0:
                    for tmp_iter in xrange(N_cuts):
                        f_traj.readline()
                        cnt_lines += 1
                        for tmp_iter_2 in xrange(Np):
                            f_index.readline()
                            f_weight.readline()
                    # else:
                    #     f_traj.readline()
                    #     cnt_lines += 1

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
                            time_past_onset_shear = cnt_lines * Delta_t_strider
                            M0 = cal_M0_simple_shear(Wi_R, box_dimension, time_past_onset_shear)
                            cnt = hist_SPE_over_beads_modified(pos, connectivity, box_dimension, hist_SPE, N_div, M0)
                            # print cnt
                        cnt_lines += 1
                    except:
                        print '[break]ing line number = ', cnt_lines
                        break
        # note that the following codes are only compatible with 3-dimenional space

    return cnt_lines


def gen_hist_SPE_arr(Np, box_dimension, N_div, R_max):
    # the range of variables are described in below:
    # r \in [0, R_max)
    # theta \in [0, \pi]
    # phi \in [0, 2\pi)
    #
    # it is of importance to follow a single definition of spherical coordinate system

    N_dimension = 3
    hist_SPE = zeros([N_div, N_div, N_div, N_dimension + 1])
    r_arr = linspace(0, R_max, N_div, endpoint=False) # endpoint=False make an open set
    dr = r_arr[1] - r_arr[0]
    
    theta_arr = linspace(0, pi, N_div) # endpoint is included because of closed set
    dtheta = theta_arr[1] - theta_arr[0]
    
    phi_arr = linspace(0, 2.*pi, N_div, endpoint=False) # endpoint=False make an open set
    dphi = phi_arr[1] - phi_arr[0]
    
    for i in range(N_div):
        r = r_arr[i] + dr/2.
        # dr/2. correct the median value of the given increment (between r_arr[i] and r_arr[i+1])
        for j in range(N_div):
            theta = theta_arr[j] + dtheta/2.
            for k in range(N_div):
                phi = phi_arr[k] + dphi/2.
                hist_SPE[i, j, k, 0] = r
                hist_SPE[i, j, k, 1] = theta
                hist_SPE[i, j, k, 2] = phi

    return hist_SPE

def get_hist_SPE_from_list(fn_list, Np, box_dimension, N_cuts, N_div, R_max, Wi_R, Delta_t_strider):
    ### This is duplicate version compared with get_hist_R_from_list
    ### Note again that the current data structure measure fdV which contains Jacobian when it has differen coordinate
    ### In consequence, the construction part can be generalized that show the same interface while contents and identification functions
    ### can be transferred through function arguments, which will be update for future.
    
    hist_SPE = gen_hist_SPE_arr(Np, box_dimension, N_div, R_max)
    
    dr = hist_SPE[1, 0, 0, 0] - hist_SPE[0, 0, 0, 0]
    dtheta = hist_SPE[0, 1, 0, 1] - hist_SPE[0, 0, 0, 1]
    dphi = hist_SPE[0, 0, 1, 2] - hist_SPE[0, 0, 1, 2]
    
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
            tmp_time_stamp_long = get_intensity_SPE_from_data(fn_base_tmp, hist_R, Np, box_dimension, N_cuts, N_div, Wi_R, Delta_t_strider)
            if (count_samples == 0):
                time_stamp_long = tmp_time_stamp_long
                Nt = time_stamp_long - N_cuts
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
        print 'USAGE of spherical coordinate:'
        print 'argv[1] == filename of list'
        print 'argv[2] == out filename'
        print 'argv[3] == number of particles'
        print 'argv[4] == box dimension'
        print 'argv[5] == N_div'
        print 'argv[6] == number of initial cuts'
        print 'argv[7] == mximally extendable length scale R_max (non-dimensional unit)'
        print 'argv[8] == Wi_R'
        print 'argv[9] == delta_t_strider = dt * strider'
    else:
        fn_list = sys.argv[1]
        fn_out = sys.argv[2]
        Np = int(sys.argv[3])
        box_dimension = float(sys.argv[4])
        # dx = float(sys.argv[5])
        N_div = float(sys.argv[5])
        N_cuts = int(sys.argv[6])
        R_max = float(sys.argv[7])
        Wi_R = float(sys.argv[8])
        Delta_t_strider = float(sys.argv[9])
        hist_R = get_hist_SPE_from_list(fn_list, Np, box_dimension, N_cuts, N_div, R_max, Wi_R, Delta_t_strider)
        save(fn_out, hist_R)
