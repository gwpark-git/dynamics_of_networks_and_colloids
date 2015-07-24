from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import mpl_toolkits.axisartist as axisartist
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from pylab import rand
import sys

from multiprocessing import Pool
from functools import partial

font_side = FontProperties()
font_side.set_size('x-small')

# from multiprocessing import Pool
# from functools import partial

def get_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
    return kd[argmin(abs(kd))] + x;

if size(sys.argv) != 9:
    print "USAGE: Plotting Trajectory Files"
    print "argv[1] == .traj"
    print "argv[2] == .hash"
    print "argv[3] == .weight"
    print "argv[4] == .ener"
    print "argv[5] == output path for figures"
    print "argv[6] == Number of particles"
    print "argv[7] == Number of skips for trajectory"
    print "argv[8] == Number of processors for parallization"
else:
    fn_traj = sys.argv[1]
    fn_connect = sys.argv[2]
    fn_weight = sys.argv[3]
    fn_ener = sys.argv[4]
    out_path = sys.argv[5]
    N_dimension = 2
    Np = int(sys.argv[6])
    N_skip = int(sys.argv[7])
    N_proc = int(sys.argv[8])

    def token(hash_data, index):
        for i in range(shape(hash_data)[1]):
            if hash_data[index][i] == "-1":
                return i
        return i

    def plot_t(given_traj, connectivity_all, color_map, dat, t):
        ft = t + 1
        # print ft
        t = t%N_proc
        # a_dat = asarray(dat)
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
        # plt.clf()
        fig = plt.figure(t, figsize=(10,6))
        # fig = plt.figure(t)
        gs = gridspec.GridSpec(2,2, width_ratios=[2, 1])
        # fig.subplots_adjust(wspace=0.4, bottom=0.3)
        ax = axisartist.Subplot(fig, gs[:,0])
        fig.add_subplot(ax)
        ax.axis([-0.2*box_dimension[0], 1.2*box_dimension[0], -0.2*box_dimension[1], 1.2*box_dimension[1]])
        ax.set_xticks(range(11))
        ax.set_yticks(range(11))
        # ax.yaxis.set_ticks_position('left')
        # ax.axis[:].major_ticks.set_tick_out(True)
        # ax.axis[:].invert_ticklabel_direction()
        # set_tick
        ax.grid('on')
        # ax.set_xlabel('x dimension')
        # ax.set_ylabel('y dimension')
        ax.set_aspect(1)
        marker_unit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]
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
        # cnt_asso = 0
        for i in range(Np):
            for j in range(Np):
                if connectivity[i,j]:
                    # cnt_asso += connectivity[i,j]
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

        annotate_text_left = 'N_dim=%ld'%(N_dimension) + '\n'
        annotate_text_left += 'N_par=%ld'%(Np) 
        ax.annotate(annotate_text_left, xy=(-0.18*box_dimension[0], 1.1*box_dimension[1]), fontsize=6, color='b', path_effects=[PathEffects.withStroke(linewidth=1, foreground='w')])

        # annotate_text_right = 'ts=%5.3e'%(given_traj[t, 0]) + '\n'
        annotate_text_right = 'ts=%ld'%((ft - 1)*N_skip) + '\n'
        annotate_text_right += 'N_bri=%ld'%(dat[ft-1, 1]) + '\n'
        annotate_text_right += 'f=%3.1f'%(dat[ft-1, 2])  
        # annotate_text_right += 'N_bri=%ld'%(cnt_asso/2) + '\n'
        # annotate_text_right += 'f=%3.1f'%(cnt_asso/Np)  

        ax.annotate(annotate_text_right, xy=(1.02*box_dimension[0], 1.1*box_dimension[1]), fontsize=6, color='b', path_effects=[PathEffects.withStroke(linewidth=1, foreground='w')])

        ax01 = axisartist.Subplot(fig, gs[0, 1])
        fig.add_subplot(ax01)
        # note that :ft on dat[:ft, 0] represent for the starting point and end point
        # ft is the real time. The reason to use ft rather than : is simple: the given dat is
        # fill with furture information since the matrix is given by reference
        # but the multiple processor is used for this computing
        mean_N_bridge = mean(dat[:ft,1])
        ref_ax01 = asarray([[dat[0, 0], mean_N_bridge],
                           [dat[ft-1, 0] + N_skip, mean_N_bridge]])
        ax01.plot(dat[:ft, 0], dat[:ft, 1], 'b.-', label = 'N_asso')
        ax01.plot(ref_ax01[:,0], ref_ax01[:,1], 'k:', linewidth=3, label = 'av.=%6.3f'%(mean_N_bridge))
        ax01.grid('on')
        ax01.axis([dat[0, 0], dat[ft-1, 0] + N_skip, 0.9*min(dat[:ft, 1]), 1.1*max(dat[:ft, 1])])
        # ax01.set_aspect(0.3)
        if ft <> 1:
            ax01.set_xticks(linspace(0, (ft - 1)*N_skip, 3))
            # range(0, ft*N_skip+2, N_skip))
        ax01.legend(loc = 'upper left', prop=font_side)
        ax01.axis["right"].toggle(ticks=False)
        # ax01.get_yaxis().tick_left()
        
        # ax01.axis[:].major_ticks.set_tick_out(True)
        # ax01.axis[:].invert_ticklabel_direction()
        # ax01.legend(loc = 'upper left')
        # ax01.plot(t, sin(t), 'b-')


        
        ax01_twin = ax01.twinx()
        ax01_twin.axis([dat[0, 0], dat[ft-1, 0] + N_skip, 0.9*min(dat[:ft, 1])*2./float(Np), 1.1*max(dat[:ft, 1])*2./float(Np)])
        
        ax11 = axisartist.Subplot(fig, gs[1, 1])
        fig.add_subplot(ax11)

        # mean_functionality = mean(dat[:ft, 2])
        mean_ener = mean(dat[:ft, 3])
        ref_ax11 = asarray([[dat[0, 0], mean_ener],
                           [dat[ft-1, 0] + N_skip, mean_ener]])
        
        ax11.plot(dat[:ft, 0], dat[:ft, 3], 'b.-', label = 'ener')
        ax11.plot(ref_ax11[:,0], ref_ax11[:,1], 'k:', linewidth=3, label = 'av.=%6.3f'%(mean_ener))
        ax11.grid('on')
        ax11.axis([dat[0, 0], dat[ft-1, 0] + N_skip, 0.9*min(dat[:ft, 3]), 1.1*max(dat[:ft, 3])])

        # ax11.set_aspect(0.3)
        # ax11.set_xticks(range(0, ft*N_skip+2, N_skip))
        if ft <> 1:
            ax11.set_xticks(linspace(0, (ft-1)*N_skip, 3))
        ax11.legend(loc = 'upper left', prop=font_side)

        # ax11.axis[:].major_ticks.set_tick_out(True)
        # ax11.axis[:].invert_ticklabel_direction()
        
        # ax11.legend(loc = 'upper left')
        # ax11.plot(t, cos(t), 'r-')
        
        # ax21 = axisartist.Subplot(fig, gs[2, 1])
        # fig.add_subplot(ax21)
        # ax21.plot([1,2,3], [4,5,6] - cos(t), 'k-')
        # ax21.set_aspect(0.3)
        # ax21.axis[:].major_ticks.set_tick_out(True)
        # ax21.axis[:].invert_ticklabel_direction()


        plt.savefig('%s/t%08d.png'%(out_path, ft-1), dpi=300, bbox_inches='tight')
        plt.close()


    if __name__ == '__main__':
        pool = Pool(processes=N_proc)
        color_map = zeros([Np, 3])
        for i in range(Np):
            color_map[i, :] = rand(3)
        cnt = 0
        dat = asarray([zeros(3)])
        with open(fn_traj, 'r') as f:
            with open(fn_connect, 'r') as f_connect:
                with open(fn_weight, 'r') as f_weight:
                    with open(fn_ener, 'r') as f_ener:
                        N_cols = 2*N_dimension*Np + 1
                        tmp_arr = zeros([N_proc, N_cols])
                        cnt_line = 0
                        c_t = arange(N_proc)
                        connectivity = zeros([N_proc, Np, Np])
                        for line in f:
                            # # st: temporal stopping
                            # if cnt > 5*N_proc:
                            #     break
                            # cnt += 1
                            # # end: temporal stopping
                            tmp_str = line.split('\t')
                            for i in range(N_cols):
                                tmp_arr[cnt_line%N_proc, i] = float(tmp_str[i])
                            hash_index = []
                            weight_index = []
                            for i in range(Np):
                                hash_index.append(f_connect.readline().split('\t')[:-1]) # last one is \n char
                                weight_index.append(f_weight.readline().split('\t')[:-1])

                            cnt_asso = 0
                            for i in range(Np):
                                index_particle = long(hash_index[i][0])
                                for j in range(1, token(hash_index, i)):
                                    index_target = long(hash_index[i][j])
                                    connectivity[cnt_line%N_proc, index_particle, index_target] = long(weight_index[i][j])
                                    cnt_asso += long(weight_index[i][j])
                            # cnt_asso = 0
                            # for i in range(Np):
                            #     for j in range(Np):
                            #         cnt_asso += connectivity[cnt_line%N_proc, i, j]
                            # tmp_dat_arr = asarray([[tmp_arr[c_t[cnt_line%N_proc], 0], cnt_asso/2., cnt_asso/float(Np)]])
                            ener = f_ener.readline().split('\t')[:-1]
                            tmp_dat_arr = asarray([[c_t[cnt_line%N_proc]*N_skip, cnt_asso/2., cnt_asso/float(Np), float(ener[1])]])

                            if cnt_line == 0:
                                dat = tmp_dat_arr
                            else:
                                dat = append(dat, tmp_dat_arr, axis=0)
                            # dat = append(dat, [[tmp_arr[c_t[cnt_line%N_proc],0], cnt_asso/2., cnt_asso/float(Np)]], axis=0)
                            # print dat
                            # print 'shape=', shape(dat)
                            # if (cnt_line <> 0 and cnt_line%N_proc == 0):
                            #     for nq in range(N_proc):
                            #         plot_t(tmp_arr, connectivity, color_map, dat, nq);
                            #     c_t += N_proc

                            cnt_line += 1
                            if (cnt_line <> 0 and cnt_line%N_proc == 0):
                                pool.map(partial(plot_t, tmp_arr, connectivity, color_map, dat), c_t)
                                c_t += N_proc
            # print dat

        # ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
        # python plot_position_multiprocessing_association_2d.py *NP40_C5900_T3.traj *NP40_C5900_T3.hash *NP40_C5900_T3.weight figures_C5900 40 12
