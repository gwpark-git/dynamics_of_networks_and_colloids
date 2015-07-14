from numpy import *
import matplotlib.pyplot as plt
import sys

N_processors = int(sys.argv[1])
index_processor = int(sys.argv[2])
N_skip = int(sys.argv[3])



traj = loadtxt('tmp.dat')
Nt = shape(traj)[0]
dimension = 2
index_x = 0
index_y = 1



Np = int((shape(traj)[1] - 1)/(2*dimension)) 
box_dimension = [1.0, 1.0]


t_block = long(Nt/N_processors)
t_st = t_block*index_processor
t_end = t_block*(index_processor + 1)

traj = traj[t_st:t_end, :]


# for t in range(t_st, t_end, N_skip):
marker_style = ['bo', 'ro', 'go', 'ko', 'co']
for t in range(0, t_block/N_skip):
    index_t = t*N_skip
    plt.clf()
    for i in range(Np):
        index_px = i*dimension*2 + 1 + index_x
        index_py = index_px + 1
        # index_vx = i*dimension*2 + 1 + dimension + index_x
        # index_vy = index_vx + 1
        plt.plot(traj[index_t, index_px], traj[index_t, index_py], marker_style[i%5])

    plt.axis('equal')
    plt.axis([0, box_dimension[0], 0, box_dimension[1]])
    plt.grid('on')
    plt.xlabel('x dimension')
    plt.ylabel('y dimension')
    plt.savefig('figures/t%08d.png'%(index_t), dpi=300)


# p = zeros([t_block, Np, dimension])
# v = zeros([t_block, Np, dimension])


# for t in range(0, t_block/N_skip):
#     index_t = t_st + t*N_skip
#     plt.clf()
#     for i in range(Np):
#         index_px = i*dimension*2 + 1 + index_x
#         index_py = index_px + 1
#         index_vx = i*dimension*2 + 1 + dimension + index_x
#         index_vy = index_vx + 1

#         for k in range(dimension):
#             # print i*dimension*2 + 1 + k, i*dimension*2 + 1 + dimension + k
#             p[t, i, k] = traj[t_st + t*N_skip, i*dimension*2 + 1 + k]%box_dimension[k]
#             v[t, i, k] = traj[t_st + t*N_skip, i*dimension*2 + 1 + dimension + k]%box_dimension[k]


# marker_style = ['bo', 'ro', 'go', 'ko', 'co']


# for t in range(t_st, t_end, N_skip):
#     plt.clf()
#     for i in range(Np):
#         plt.plot(p[t][i,0], p[t][i,1], marker_style[i%5])

    # plt.savefig('figures/t%08d.pdf'%(t))

