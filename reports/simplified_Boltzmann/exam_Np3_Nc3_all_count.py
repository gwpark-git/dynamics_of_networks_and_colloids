

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
        # try:
        #     if not dat[i, 1] in p:
        #         # print dat[i,1]
        #         # cat.append(dat[i,:])
        #         ind.append(int(dat[i, 0]))
        #         p.append(dat[i, 1])
        # except:
        #     ind.append(int(dat[i, 0]))
        #     p.append(dat[i, 1])
            # cat.append(dat[i,:])
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

# r0 = asarray([norm(p0 - p0), norm(p0 - p1), norm(p0 - p2)])
# r1 = asarray([norm(p1 - p1), norm(p1 - p2), norm(p1 - p0)])
# r2 = 

Pr = exp(-r)

Np = 3
Nc = 3
# ((Np, Nc)) = (Np + Nc - 1, Nc) => 27!/(25!2!) = 27*13 = 351
# Ns = 21

# w0 = zeros([Ns, Np])
# w1 = zeros([Ns, Np])
# w2 = zeros([Ns, Np])

# def get_combination_array(Np, Nc):
#     Ns = comb(Np + Nc - 1, Nc)
#     # w = zeros([Ns, Np])
#     for s in range(Ns):
#         for i in range(Np):

# w = asarray([[5, 0, 0],
#              [4, 1, 0], [4, 0, 1],
#              [3, 2, 0], [3, 1, 1], [3, 0, 2],
#              [2, 3, 0], [2, 2, 1], [2, 1, 2], [2, 0, 3],
#              [1, 4, 0], [1, 3, 1], [1, 2, 2], [1, 1, 3], [1, 0, 4],
#              [0, 5, 0], [0, 4, 1], [0, 3, 2], [0, 2, 3], [0, 1, 4], [0, 0, 5]])
Ns_distinguishable = 4
Ns_indist = 6
u_indist = asarray([uij(norm(p0-p0)), uij(norm(p1-p1)), uij(norm(p2-p2)), uij(norm(p0-p1)), uij(norm(p0-p2)), uij(norm(p1-p2))])
# degeneracy = asarray([3, 1, 1, 1])
degeneracy_roof = Np
Nc_association = Np*Nc
w = get_weight_combinations(Ns_indist, Nc_association)
Ns = shape(w)[0]
pot = zeros(Ns)
for s in range(Ns):
    for i in range(Ns_indist):
        pot[s] += w[s, i]*u_indist[i]
tmp_P = exp(-pot)
Z = sum(tmp_P)
P = tmp_P/Z

lower_lim = 1e-6
dat = []
for i,p in enumerate(P):
    if p > lower_lim:
        dat.append([i,p])
dat = asarray(dat)


diff_P_ind, diff_P_probability = find_diff_P(dat, lower_lim/2.0)
diff_levels = size(diff_P_ind)
cnt_levels = zeros(diff_levels)
for i in range(diff_levels):
    # ref_P = diff_P_probability[i]
    for j in range(shape(dat)[0]):
        if abs(dat[j, 1] - diff_P_probability[i]) < lower_lim/2.0:
            cnt_levels[i] += 1

lim_x = 1.5*Ns
ref_P = asarray([[0, 1],
                 [Ns, 1]])

fontP = FontProperties(size='small')
import matplotlib.colors as colors
import matplotlib.cm as cmx
cm = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=diff_levels)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

plt.clf()
plt.ion()
xmin, xmax, ymin, ymax = 0, 1.5*Ns, lower_lim, 10**-1
# D_y = (ymax/ymin)/100000.
D_y = 1.1
shift_y = 1.2
rep_p = 0.
av_NAS = 0.
for i in range(diff_levels):
    rep_p += cnt_levels[i]*diff_P_probability[i]
    av_NAS += cnt_levels[i]*diff_P_probability[i]*sum(w[diff_P_ind[i]][degeneracy_roof:])
plt.semilogy(dat[:,0], dat[:,1], 'b.', label = 'all the state\n(u11 + u22 + u33, u12, u13, u23)')
for i in range(diff_levels):
    color_val = scalarMap.to_rgba(range(diff_levels)[i])
    plt.semilogy(ref_P[:,0], ref_P[:,1]*diff_P_probability[i], '-', color=color_val, label = '%d level: %4.3e\n%ld times of (%d, %d, %d, %d)\noccupied_P = %4.6lf'%(i, diff_P_probability[i], cnt_levels[i], w[diff_P_ind[i]][0] + w[diff_P_ind[i]][1] + w[diff_P_ind[i]][2], w[diff_P_ind[i]][3], w[diff_P_ind[i]][4], w[diff_P_ind[i]][5], cnt_levels[i]*diff_P_probability[i]))
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
    plt.annotate('%d level'%(i), xy=(Ns, y_pos))
plt.legend(loc = 'upper right', prop=fontP, numpoints=1)
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel('index for state')
plt.ylabel('log10 scale for probability of specific state')
plt.title('expected NAS = %4.3lf out of %d chains (%4.6lf ratio), rep_p = %4.6lf'%(av_NAS, Nc*Np, av_NAS/(Nc*Np), rep_p))
plt.show()

    # return asarray(cat)

# tmp_P0 = zeros(Ns)
# tmp_P1 = zeros(Ns)
# tmp_P2 = zeros(Ns)
# for i in range(Ns):
#     for j in range(Np):
#         tmp_P0[i] += w[i,j]*uij(r[0, j])
#         tmp_P1[i] += w[i,j]*uij(r[1, j])
#         tmp_P2[i] += w[i,j]*uij(r[2, j])
#     tmp_P0[i] = exp(-tmp_P0[i])
#     tmp_P1[i] = exp(-tmp_P1[i])
#     tmp_P2[i] = exp(-tmp_P2[i])

# Z0, Z1, Z2 = sum(tmp_P0), sum(tmp_P1), sum(tmp_P2)
# P0, P1, P2 = tmp_P0/Z0,   tmp_P1/Z1,   tmp_P2/Z2

# plt.clf()
# plt.ion()
# plt.plot(P0, 'b.-', label='P0, max=(%d, %4.3f)'%(argmax(P0), P0[argmax(P0)]))
# plt.plot(P1, 'r.-', label='P1, max=(%d, %4.3f)'%(argmax(P1), P1[argmax(P1)]))
# plt.plot(P2, 'g.-', label='P2, max=(%d, %4.3f)'%(argmax(P2), P2[argmax(P2)]))
# plt.grid()
# plt.legend(loc = 'upper left')
# plt.xlabel('index for microstate')
# plt.ylabel('probability for microstate')
# plt.axis([0, Ns-1, -0.01, 1.0])
# plt.show()
