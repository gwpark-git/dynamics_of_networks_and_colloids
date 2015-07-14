from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from multiprocessing import Pool
from functools import partial

fn = sys.argv[1]
out_path = sys.argv[2]

def plot_t(traj, t):
    dimension = 3
    Np = int((shape(traj)[1] - 1)/(2*dimension))
    box_dimension = [10.0, 10.0, 10.0]

    color_style = ['b', 'r', 'g', 'k', 'c', 'm', 'purple', 'cyan' ]
    marker_style = '.'
    fig = plt.figure(t)
    ax = fig.add_subplot(111, projection='3d')
    ax.grid('on')
    ax.set_xlabel('x dimension')
    ax.set_ylabel('y dimension')
    ax.set_zlabel('z dimension')
    ax.set_aspect(1)
    # marker_unit = (ax.transData.transform((1,0,0)) - ax.transData.transform((0,0,0)))[0]
    for i in range(Np):
        index_px = i*dimension*2 + 1 + 0
        index_py = index_px + 1
        index_pz = index_py + 1
        ax.scatter(traj[t, index_px], traj[t, index_py], traj[t, index_pz], color = color_style[i%5], marker=marker_style)

    # ax.axis('equal')
    # ax.axis([0, box_dimension[0], 0, box_dimension[1], 0, box_dimension[2]])
    # plt.axis([-0.5*box_dimension[0], 0.5*box_dimension[0], -0.5*box_dimension[1], 0.5*box_dimension[1]])
    plt.savefig('%s/t%08d.png'%(out_path, t), dpi=300, bbox_inches = 'tight')
    plt.close()


if __name__ == '__main__':
    traj = loadtxt(fn)
    # Nt = 1000000
    Nt = shape(traj)[0]
    N_skip = 1 #already skipped data
    pool = Pool(processes=None) # processes = None allocate each thread to the each CPUs.
    pool.map(partial(plot_t, traj), range(0, Nt, N_skip))

# ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
