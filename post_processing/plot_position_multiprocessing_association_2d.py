from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from pylab import rand
import sys

from multiprocessing import Pool
from functools import partial

# from multiprocessing import Pool
# from functools import partial

def get_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
    return kd[argmin(abs(kd))] + x;

if size(sys.argv) != 7:
    print "USAGE: Plotting Trajectory Files"
    print "argv[1] == .traj"
    print "argv[2] == .hash"
    print "argv[3] == .weight"
    print "argv[4] == output path for figures"
    print "argv[5] == Number of particles"
    print "argv[6] == Number of processors for parallization"
else:
    fn_traj = sys.argv[1]
    fn_connect = sys.argv[2]
    fn_weight = sys.argv[3]
    out_path = sys.argv[4]
    N_dimension = 2
    Np = int(sys.argv[5])
    N_proc = int(sys.argv[6])


    # given_traj = loadtxt(fn)
    # hash_index = loadtxt('NP80_C100_T3.hash')
    # weight_index = loadtxt('NP80_C100_T3.weight')


    # connectivity = zeros([Np, Np])
    # weight = zeros([Np, Np])

    def token(hash_data, index):
        for i in range(shape(hash_data)[1]):
            if hash_data[index][i] == "-1":
                return i
        return i

    # for i in range(Np):
    #     index_particle = hash_index[i, 0]
    #     for j in range(1, token(hash_index, i)):
    #         index_target = hash_index[i, j]
    #         connectivity[index_particle, index_target] = weight_index[i,j]

    def plot_t(given_traj, connectivity_all, color_map, t):
        ft = t
        t = t%N_proc
        connectivity = connectivity_all[t, :, :]
        box_dimension = [10.0, 10.0]

        ref_PBC_left = asarray([[0, 0],
                                [0, 10]])
        ref_PBC_right = asarray([[10, 0],
                                 [10, 10]])
        ref_PBC_bottom = asarray([[0, 0],
                                  [10, 0]])
        ref_PBC_top = asarray([[0, 10],
                               [10, 10]])

        marker_style = 'o'
        plt.clf()
        fig = plt.figure(t)
        ax = fig.add_subplot(111)
        ax.axis([-0.2*box_dimension[0], 1.2*box_dimension[0], -0.2*box_dimension[1], 1.2*box_dimension[1]])
        ax.set_xticks(range(11))
        ax.set_yticks(range(11))
        ax.grid('on')
        ax.set_xlabel('x dimension')
        ax.set_ylabel('y dimension')
        ax.set_aspect(1)
        marker_unit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]/1.4
        for i in range(Np):
            index_px = i*N_dimension*2 + 1 + 0
            index_py = index_px + 1
            ax.plot(given_traj[t, index_px], given_traj[t, index_py], color=color_map[i,:], marker=marker_style, markersize=marker_unit, alpha=0.8)

        # note that following code is only valid for 2-dimensional space
        # for general purpose, it should be changed by recursive call
            for shift_fac_x in [-1., 0., 1.]:
                for shift_fac_y in [-1., 0., 1.]:
                    if not(shift_fac_x == 0. and shift_fac_y == 0.):
                        ax.plot(given_traj[t, index_px] + shift_fac_x*box_dimension[0], given_traj[t, index_py] + shift_fac_y*box_dimension[1], color=color_map[i,:], marker=marker_style, markersize=marker_unit, alpha=0.4)


        tmp_arr = zeros([2, N_dimension])
        cnt_asso = 0
        for i in range(Np):
            for j in range(Np):
                if connectivity[i,j]:
                    cnt_asso += connectivity[i,j]
                    for k in range(N_dimension):
                        index_pi_k = i*N_dimension*2 + 1 + k
                        index_pj_k = j*N_dimension*2 + 1 + k
                        tmp_arr[0, k] = given_traj[t, index_pi_k]
                        tmp_arr[1, k] = get_minimum_distance_k_from_x(given_traj[t, index_pi_k], given_traj[t, index_pj_k], box_dimension[k])
                    # tmp_arr[1, k] = given_traj[t, index_pj_k]
                    ax.plot(tmp_arr[:,0], tmp_arr[:,1], 'r-', linewidth=0.5)

                    mean_x = 0.5*(tmp_arr[0, 0] + tmp_arr[1, 0])
                    mean_y = 0.5*(tmp_arr[0, 1] + tmp_arr[1, 1])
                    ax.plot(tmp_arr[:,0], tmp_arr[:,1], 'k.', markersize=3)
                    if connectivity[i,j] > 1:
                        ax.annotate('%d'%(connectivity[i, j]), xy=(mean_x, mean_y), fontsize=4, color='b', path_effects=[PathEffects.withStroke(linewidth=1, foreground='w')])


        ax.plot(ref_PBC_left[:,0], ref_PBC_left[:,1], 'k--', linewidth=3)
        ax.plot(ref_PBC_right[:,0], ref_PBC_right[:,1], 'k--', linewidth=3)
        ax.plot(ref_PBC_bottom[:,0], ref_PBC_bottom[:,1], 'k--', linewidth=3)
        ax.plot(ref_PBC_top[:,0], ref_PBC_top[:,1], 'k--', linewidth=3)

        plt.title('(N_association, N_beads, f) = (%ld, %ld, %6.1f)'%(cnt_asso/2, Np, cnt_asso/(2.*Np)))
        plt.savefig('%s/t%08d.png'%(out_path, ft), dpi=300, bbox_inches='tight')
        plt.close()
    # plt.show()
    # plt.savefig('%s/t%08d.png'%(out_path, ft), dpi=300, bbox_inches = 'tight')
    # plt.close()
    # plt.show()


    if __name__ == '__main__':
        pool = Pool(processes=N_proc)
        color_map = zeros([Np, 3])
        for i in range(Np):
            color_map[i, :] = rand(3)

        with open(fn_traj, 'r') as f:
            with open(fn_connect, 'r') as f_connect:
                with open(fn_weight, 'r') as f_weight:
                    N_cols = 2*N_dimension*Np + 1
                    tmp_arr = zeros([N_proc, N_cols])
                    cnt_line = 0
                    c_t = arange(N_proc)
                    connectivity = zeros([N_proc, Np, Np])
                    for line in f:
                        tmp_str = line.split('\t')
                        for i in range(N_cols):
                            tmp_arr[cnt_line%N_proc, i] = float(tmp_str[i])
                        hash_index = []
                        weight_index = []
                        for i in range(Np):
                            hash_index.append(f_connect.readline().split('\t')[:-1]) # last one is \n char
                            weight_index.append(f_weight.readline().split('\t')[:-1])

                        for i in range(Np):
                            index_particle = long(hash_index[i][0])
                            for j in range(1, token(hash_index, i)):
                                index_target = long(hash_index[i][j])
                                connectivity[cnt_line%N_proc, index_particle, index_target] = long(weight_index[i][j])


                        cnt_line += 1
                        if (cnt_line <> 0 and cnt_line%N_proc == 0):
                            pool.map(partial(plot_t, tmp_arr, connectivity, color_map), c_t)
                            c_t += N_proc


    # ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
    # python plot_position_multiprocessing_association_2d.py *NP40_C5900_T3.traj *NP40_C5900_T3.hash *NP40_C5900_T3.weight figures_C5900 40 12
