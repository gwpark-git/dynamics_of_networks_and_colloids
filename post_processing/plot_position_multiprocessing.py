from numpy import *
import matplotlib.pyplot as plt
import sys
from multiprocessing import Pool
from functools import partial

fn = sys.argv[1]
out_path = sys.argv[2]

def plot_t(traj, t):
    dimension = 2
    Np = int((shape(traj)[1] - 1)/(2*dimension))
    L_b = 10.0
    box_dimension = [L_b, L_b]

    marker_style = ['bo', 'ro', 'go', 'ko', 'co']
    fig = plt.figure(t)
    ax = fig.add_subplot(111)
    ax.axis([0, box_dimension[0], 0, box_dimension[1]])
    ax.grid('on')
    ax.set_xlabel('x dimension')
    ax.set_ylabel('y dimension')
    ax.set_aspect(1)
    # print ax.transData.transform((0,1))
    # print ax.transData.transform((0,0))
    # the following transfer size of marker as 1. That value is related with the pixel between 0 and 1.
    # marker_unit, = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0))
    tmp_dat, marker_unit = (ax.transData.transform((0, 1)) - ax.transData.transform((0,0)))


    for i in range(Np):
        index_px = i*dimension*2 + 1 + 0
        index_py = index_px + 1
        ax.plot(traj[t, index_px], traj[t, index_py], marker_style[i%5], markersize=marker_unit)

    # ax.axis('equal')
    # plt.axis([-0.5*box_dimension[0], 0.5*box_dimension[0], -0.5*box_dimension[1], 0.5*box_dimension[1]])
    plt.savefig('%s/t%08d.png'%(out_path, t), dpi=300, bbox_inches = 'tight')
    plt.close()


if __name__ == '__main__':
    traj = loadtxt(fn)
    # Nt = 1000000
    Nt = shape(traj)[0]
    N_skip = 1
    pool = Pool(processes=None) # processes = None allocate each thread to the each CPUs.
    pool.map(partial(plot_t, traj), range(0, Nt, N_skip))

# ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
