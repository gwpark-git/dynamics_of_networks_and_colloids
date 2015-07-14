from numpy import *
import matplotlib.pyplot as plt
import sys
from multiprocessing import Pool
from functools import partial

def plot_t(traj, t):
    dimension = 2
    Np = int((shape(traj)[1] - 1)/(2*dimension)) 
    box_dimension = [10.0, 10.0]

    marker_style = ['bo', 'ro', 'go', 'ko', 'co']
    fig = plt.figure(t)
    ax = fig.add_subplot(111)
    for i in range(Np):
        index_px = i*dimension*2 + 1 + 0
        index_py = index_px + 1
        ax.plot(traj[t, index_px], traj[t, index_py], marker_style[i%5])

    # ax.axis('equal')
    ax.axis([0, box_dimension[0], 0, box_dimension[1]])
    # plt.axis([-0.5*box_dimension[0], 0.5*box_dimension[0], -0.5*box_dimension[1], 0.5*box_dimension[1]])
    ax.grid('on')
    ax.set_xlabel('x dimension')
    ax.set_ylabel('y dimension')
    ax.set_aspect(1)
    plt.savefig('test_t%08d.png'%(t), dpi=300)
    plt.close()


if __name__ == '__main__':
    traj = loadtxt('tmp.dat')
    # Nt = 1000000
    Nt = shape(traj)[0]
    N_skip = 1 #already skipped data
    plot_t(traj, 0)
    # pool = Pool(processes=None) # processes = None allocate each thread to the each CPUs.
    # pool.map(partial(plot_t, traj), range(0, Nt, N_skip))
