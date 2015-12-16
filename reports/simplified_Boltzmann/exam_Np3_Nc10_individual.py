

from numpy import *
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.linalg import norm
from scipy.misc import comb
from combinatorial_sequence import *

def find_diff_P(dat, tolerance):
    ind = []
    p = []
    Nd = shape(dat)[0]
    # Nd = 5
    for i in range(Nd):
        ident = 1
        for ref in p:
            if abs(dat[i, 1] - ref) < tolerance:
                ident = 0
        if ident:
            ind.append(dat[i, 0])
            p.append(dat[i, 1])
    sorted_index = argsort(p)
    re_ind = zeros(size(ind))
    re_p = zeros(size(ind))
    for i in range(size(ind)):
        re_ind[i] = ind[sorted_index[i]]
        re_p[i] = p[sorted_index[i]]
    return re_ind[::-1], re_p[::-1]


def uij(distance):
    return distance**2.0

p0 = asarray([0, 0])
p1 = asarray([1, 1])
p2 = asarray([0, -2])

r = asarray([[0, norm(p0 - p1), norm(p0 - p2)],
             [norm(p1 - p0), 0, norm(p1 - p2)],
             [norm(p2 - p0), norm(p2 - p1), 0]])


Pr = exp(-r)

Np = 3
Nc = 25
Ns = asarray([3, 3, 3])
Ns_indist = 3
u_indist = asarray([[uij(norm(p0-p0)), uij(norm(p0-p1)), uij(norm(p0-p2))],
                    [uij(norm(p1-p1)), uij(norm(p1-p2)), uij(norm(p1-p0))],
                    [uij(norm(p2-p2)), uij(norm(p2-p0)), uij(norm(p2-p1))]])
w = get_weight_combinations(Np, Nc)
# w = asarray([get_weight_combinations(Ns[0], Nc), get_weight_combinations(Ns[1], Nc), get_weight_combination(Ns[2], Nc)])
# Ns_distinguishable = 4
# Ns_indist = 6
# u_indist = asarray([uij(norm(p0-p0)), uij(norm(p1-p1)), uij(norm(p2-p2)), uij(norm(p0-p1)), uij(norm(p0-p2)), uij(norm(p1-p2))])
# degeneracy_roof = Np
# Nc_association = Np*Nc
# w = get_weight_combinations(Ns_indist, Nc_association)
Ns = shape(w)[0]
pot = zeros([Ns, Np])
for k in range(Np):
    for s in range(Ns):
        for i in range(Ns_indist):
            pot[s, k] += w[s, i]*u_indist[k,i]
tmp_P = exp(-pot)
Z = zeros(Np)
P = zeros([Ns, Np])
for i in range(Np):
    Z[i] = sum(tmp_P[:,i])
    P[:, i] = tmp_P[:,i]/Z[i]

P_map = ones(Ns**Np)
P_index = zeros([Ns**Np, 3])
for i in range(Ns):
    for j in range(Ns):
        for k in range(Ns):
            index_now = k + j*Ns + i*Ns*Ns
            P_map[index_now] = P[i,0]*P[j,1]*P[k, 2]
            P_index[index_now, :] = asarray([i, j, k])

lower_lim = 1e-6

dat = []
for i, p in enumerate(P_map):
    if p > lower_lim:
        dat.append([i, p])
dat = asarray(dat)

diff_P_ind, diff_P_probability = find_diff_P(dat, lower_lim/2.0)
diff_levels = size(diff_P_ind)
cnt_levels = zeros(diff_levels)
for i in range(diff_levels):
    # ref_P = diff_P_probability[i]
    for j in range(shape(dat)[0]):
        if abs(dat[j, 1] - diff_P_probability[i]) < lower_lim/2.0:
            cnt_levels[i] += 1

rep_p = 0.
av_NAS = 0.
for i in range(diff_levels):
    rep_p += cnt_levels[i]*diff_P_probability[i]
    for j in range(Np):
        av_NAS += cnt_levels[i]*diff_P_probability[i]*sum(w[P_index[diff_P_ind[i],j]][1:])

lim_x = 1.5*Ns
ref_P = asarray([[0, 1],
                 [Ns**Np, 1]])



fontP = FontProperties(size='x-small')
import matplotlib.colors as colors
import matplotlib.cm as cmx
cm = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=diff_levels)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

plt.clf()
plt.ion()
xmin, xmax, ymin, ymax = 0, 1.8*Ns**Np, lower_lim, 10**0
D_y = 1.1
shift_y = 1.2
plt.semilogy(dat[:,0], dat[:,1], 'b.', label = 'all the state\n(u11+u22+u33,u12,u13,u23)')

for i in range(diff_levels):
    color_val = scalarMap.to_rgba(range(diff_levels)[i])
    plt.semilogy(ref_P[:,0], ref_P[:,1]*diff_P_probability[i], '-', color=color_val, label = '%d level: %4.3e\n%ld times of (%d, %d, %d, %d)\noccupied_P = %4.6lf'%(i, diff_P_probability[i], cnt_levels[i], w[P_index[diff_P_ind[i],0]][0] + w[P_index[diff_P_ind[i],1]][0] + w[P_index[diff_P_ind[i],2]][0], w[P_index[diff_P_ind[i],0]][1] + w[P_index[diff_P_ind[i],1]][2], w[P_index[diff_P_ind[i],0]][2] + w[P_index[diff_P_ind[i], 2]][1], w[P_index[diff_P_ind[i], 1]][1] + w[P_index[diff_P_ind[i], 2]][2], cnt_levels[i]*diff_P_probability[i]))
    p=0
    if i==0:
        y_pos = diff_P_probability[i]
        p=0
    elif y_pos/diff_P_probability[i] < D_y:
        p += 1        
        y_pos = diff_P_probability[i]/(p*shift_y)
    else:
        y_pos = diff_P_probability[i]
        p=0
    plt.annotate('%d level'%(i), xy=(Ns**Np, y_pos))
plt.legend(loc = 'upper right', prop=fontP, numpoints=1)
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel('index for state')
plt.ylabel('log10 scale for probability of specific state')
plt.title('expected NAS = %4.3lf out of %d chains (%4.6lf ratio), rep_p = %4.6lf'%(av_NAS, Nc*Np, av_NAS/(Nc*Np), rep_p))
plt.show()
