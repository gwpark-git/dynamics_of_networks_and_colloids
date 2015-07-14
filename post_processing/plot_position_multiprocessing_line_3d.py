from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from multiprocessing import Pool
from functools import partial

fn = sys.argv[1]
out_path = sys.argv[2]
N_proc = int(sys.argv[3])
N_dimension = 3
Np = 1000


def plot_t(given_traj, t):
    ft = t
    t = t%N_proc
    
    # Np = int((shape(traj)[1] - 1)/(2*dimension)) 
    box_dimension = [10.0, 10.0, 10.0]

    color_style = ['b', 'r', 'g', 'k', 'c', 'm', 'purple', 'cyan' ]
    marker_style = '.'
    fig = plt.figure(t)
    ax = fig.add_subplot(111, projection='3d')
    for i in range(Np):
        index_px = i*N_dimension*2 + 1 + 0
        index_py = index_px + 1
        index_pz = index_py + 1
        ax.scatter(given_traj[t, index_px], given_traj[t, index_py], given_traj[t, index_pz], color = color_style[i%5], marker=marker_style)

    ax.grid('on')
    ax.set_xlabel('x dimension')
    ax.set_ylabel('y dimension')
    ax.set_zlabel('z dimension')
    ax.set_aspect(1)
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
            if (cnt_line <> 0 and cnt_line%N_proc == 0):
                pool.map(partial(plot_t, tmp_arr), c_t)
                c_t += N_proc

            cnt_line += 1
        
# ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
