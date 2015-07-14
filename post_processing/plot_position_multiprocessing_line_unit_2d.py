from numpy import *
import matplotlib.pyplot as plt
import sys
from multiprocessing import Pool
from functools import partial

fn = sys.argv[1]
out_path = sys.argv[2]
N_dimension = 2
Np = int(sys.argv[3])
N_proc = int(sys.argv[4])


def plot_t(given_traj, t):
    ft = t
    t = t%N_proc
    box_dimension = [10.0, 10.0]

    color_style = ['b', 'r', 'g', 'k', 'c', 'm', 'purple', 'cyan' ]
    marker_style = 'o'
    fig = plt.figure(t)
    ax = fig.add_subplot(111)
    ax.axis([0, box_dimension[0], 0, box_dimension[1]])
    ax.set_xticks(range(11))
    ax.set_yticks(range(11))
    ax.grid('on')
    ax.set_xlabel('x dimension')
    ax.set_ylabel('y dimension')
    ax.set_aspect(1)
    marker_unit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]
    for i in range(Np):
        index_px = i*N_dimension*2 + 1 + 0
        index_py = index_px + 1
        ax.plot(given_traj[t, index_px], given_traj[t, index_py], color=color_style[i%size(color_style)], marker=marker_style, markersize=marker_unit)
        
    plt.savefig('%s/t%08d.png'%(out_path, ft), dpi=300, bbox_inches = 'tight')
    plt.close()


if __name__ == '__main__':
    pool = Pool(processes=N_proc)
    with open(fn, 'r') as f:
        N_cols = 2*N_dimension*Np + 1
        tmp_arr = zeros([N_proc, N_cols])
        cnt_line = 0
        c_t = arange(N_proc)
        for line in f:
            tmp_str = line.split('\t')
            for i in range(N_cols):
                tmp_arr[cnt_line%N_proc, i] = float(tmp_str[i])
            cnt_line += 1
            if (cnt_line <> 0 and cnt_line%N_proc == 0):
                pool.map(partial(plot_t, tmp_arr), c_t)
                c_t += N_proc
        

# ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
