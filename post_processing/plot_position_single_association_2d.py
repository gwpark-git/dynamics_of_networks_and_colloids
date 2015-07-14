from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from pylab import rand
import sys
# from multiprocessing import Pool
# from functools import partial

def get_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
    return kd[argmin(abs(kd))] + x;

fn = 'NP40_C100_T3.traj'
out_path = ''
N_dimension = 2
Np = 40
# N_proc = int(sys.argv[4])


# def plot_t(given_traj, t):
#     ft = t
#     t = t%N_proc
given_traj = loadtxt(fn)
hash_index = loadtxt('NP40_C100_T3.hash')
weight_index = loadtxt('NP40_C100_T3.weight')


connectivity = zeros([Np, Np])
# weight = zeros([Np, Np])

def token(hash_data, index):
    for i in range(shape(hash_data)[1]):
        if hash_data[index, i] == -1:
            return i
    return i

for i in range(Np):
    index_particle = hash_index[i, 0]
    for j in range(1, token(hash_index, i)):
        index_target = hash_index[i, j]
        connectivity[index_particle, index_target] = weight_index[i,j]

t = 0
ft = t
box_dimension = [10.0, 10.0]

ref_PBC_left = asarray([[0, 0],
                        [0, 10]])
ref_PBC_right = asarray([[10, 0],
                         [10, 10]])
ref_PBC_bottom = asarray([[0, 0],
                          [10, 0]])
ref_PBC_top = asarray([[0, 10],
                       [10, 10]])

# color_style = ['y', 'g', 'c', 'm', 'cyan', ]
# color_style = ['g', 'y', 'cyan']
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
color_map = zeros([Np, 3])
for i in range(Np):
    color_map[i, :] = rand(3)
marker_unit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]/1.3
# marker_unit = diff(ax.transData.transform(zip([0]*len(dx), dx)))
# marker_unit = ax.transData.transform([0])
for i in range(Np):
    index_px = i*N_dimension*2 + 1 + 0
    index_py = index_px + 1
    # plt.Circle(given_traj[t, index_px], given_traj[t, index_py], color=color_style[i%size(color_style)], radius=None, fill=False)
    ax.plot(given_traj[t, index_px], given_traj[t, index_py], color=color_map[i,:], marker=marker_style, markersize=marker_unit)
    
    # note that following code is only valid for 2-dimensional space
    # for general purpose, it should be changed by recursive call
    for shift_fac_x in [-1., 0., 1.]:
        for shift_fac_y in [-1., 0., 1.]:
            if not(shift_fac_x == 0. and shift_fac_y == 0.):
                ax.plot(given_traj[t, index_px] + shift_fac_x*box_dimension[0], given_traj[t, index_py] + shift_fac_y*box_dimension[1], color=color_map[i,:], marker=marker_style, markersize=marker_unit)
            

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
            ax.plot(tmp_arr[:,0], tmp_arr[:,1], 'r.-', linewidth=0.5, markersize=3)
            mean_x = 0.5*(tmp_arr[0, 0] + tmp_arr[1, 0])
            mean_y = 0.5*(tmp_arr[0, 1] + tmp_arr[1, 1])
            ax.annotate('%d'%(connectivity[i, j]), xy=(mean_x, mean_y), fontsize=6, color='b', path_effects=[PathEffects.withStroke(linewidth=2, foreground='w')])


ax.plot(ref_PBC_left[:,0], ref_PBC_left[:,1], 'k--', linewidth=3)
ax.plot(ref_PBC_right[:,0], ref_PBC_right[:,1], 'k--', linewidth=3)
ax.plot(ref_PBC_bottom[:,0], ref_PBC_bottom[:,1], 'k--', linewidth=3)
ax.plot(ref_PBC_top[:,0], ref_PBC_top[:,1], 'k--', linewidth=3)

plt.title('N_association = %ld'%(cnt_asso/2))
plt.savefig('tmp.pdf', dpi=300, bbox_inches='tight')
# plt.show()
# plt.savefig('%s/t%08d.png'%(out_path, ft), dpi=300, bbox_inches = 'tight')
# plt.close()
# plt.show()

# if __name__ == '__main__':
#     pool = Pool(processes=N_proc)
# with open(fn, 'r') as f:
#     N_cols = 2*N_dimension*Np + 1
#     tmp_arr = zeros([N_proc, N_cols])
#     cnt_line = 0
#     c_t = arange(N_proc)
#     for line in f:
#         tmp_str = line.split('\t')
#         for i in range(N_cols):
#             tmp_arr[cnt_line%N_proc, i] = float(tmp_str[i])
#         cnt_line += 1
#         if (cnt_line <> 0 and cnt_line%N_proc == 0):
#             plot_t(tmp_arr, c_t)
#             # pool.map(partial(plot_t, tmp_arr), c_t)
#             c_t += N_proc
        

# ffmpeg -r 60 -i figures/t%08d.png -vcodec copy out.mov
