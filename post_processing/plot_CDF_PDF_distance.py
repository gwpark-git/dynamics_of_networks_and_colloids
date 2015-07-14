

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
import sys
from chebyshev import *
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('x-small')

# def inv_map_diff_from_polynomial(x, c):
#     dG_F = diff_gen_from_polynomial(x, c)
#     dF_r = 1./dG_F
#     return dF_r

# x = loadtxt('cdf_r.dat')
# y = arange(0, size(x),1)/float(size(x))
x = sort(loadtxt(sys.argv[1])[:,5])
# x = sort(loadtxt('100.info')[:,5])
cut_min = float(sys.argv[4])
cut_max = float(sys.argv[5])
for cnt in range(size(x)):
    if x[cnt] > cut_min:
        break;
x = x[cnt:]
    
for cnt in range(size(x)):
    if x[cnt] > cut_max:
        break;
x = x[:cnt]

t = []
for i in linspace(0,size(x)-1, 5000):
    t.append(x[int(i)])
N_t = size(t)
y = arange(0, N_t, 1)/float(N_t)
t_min = min(t)
t_max = max(t)
tck = interpolate.splrep(t,y,k=3,s=0)
tnew = linspace(t_min, t_max, 5000)
# xnew = arange(min(x), max(x), 0.0001)
ynew = interpolate.splev(tnew, tck, der=0)
diff = interpolate.splev(tnew, tck, der=1)



# t_lim = tnew
# y_lim = ynew
c = get_polynomial_from_Cheby(tnew, ynew, int(sys.argv[3]))
t_cheby = linspace(min(tnew), max(tnew), 100)
y_cheby = gen_from_polynomial(t_cheby, c)
diff_cheby = diff_gen_from_polynomial(t_cheby, c)

# max_index_diff_cheby = argmax(diff_cheby)
ref_max = asarray([[t_cheby[argmax(diff_cheby)], -0.1],
                   [t_cheby[argmax(diff_cheby)], max(diff_cheby)+0.1]])


ref_line = asarray([[1., -0.1],
                    [1., max(diff_cheby)+0.1]])

# max_diff_cheby = max(diff_cheby)
# N_diff = size(diff_cheby)
# for i in range(N_diff):
#     if N_diff

# cnt = 0
# for cnt in range(size(tnew)):
#     if tnew[cnt] > 1.:
#         break;
# t_lim = t
# y_lim = y
# c = get_polynomial_from_Cheby(y_lim, t_lim, 12)
# y_cheby = linspace(min(y_lim), max(y_lim), 100)
# t_cheby = gen_from_polynomial(y_cheby, c)
# # t_cheby = linspace(min(t_lim), max(t_lim), 100)
# # y_cheby = gen_from_polynomial(t_cheby, c)
# # diff_cheby = diff_gen_from_polynomial(y_cheby, c)
# diff_cheby = inv_map_diff_from_polynomial(y_cheby, c)


plt.clf()
fig, ax1 = plt.subplots()
# leg = []
ax2 = ax1.twinx()
leg_origin = []

leg_origin.append(ax1.plot(t, y, 'bo', markerfacecolor='white', markeredgecolor='blue', label='CDF(r), data')[0])

leg_origin.append(ax1.plot(tnew, ynew, 'g-', linewidth=2, label='CDF(r), pointwise')[0])
leg_origin.append(ax2.plot(tnew, diff, 'g-', label = 'PDF(r), pointwise')[0])
leg_origin.append(ax2.plot(ref_line[:,0], ref_line[:,1], 'k--', linewidth=3, label='Range for Repulsion')[0])

leg_origin.append(ax2.plot(ref_max[:,0], ref_max[:,1], 'r--', linewidth=3, label = 'max_pdf_r =  %6.3e'%(t_cheby[argmax(diff_cheby)]))[0])

leg_origin.append(ax1.plot(t_cheby, y_cheby, 'r.', label='CDF(r), Chebyshev')[0])
# leg.append(ax2.plot(tnew, 4.*exp(-U_ij(tnew)), 'g-', label = 'U(r)')[0])

leg_origin.append(ax2.plot(t_cheby, diff_cheby, 'r-', linewidth=3, label='PDF(r), Chebyshev')[0])

ax1.grid('on')
ax1.set_xlabel('distance, r')
ax1.set_ylabel('cumulative distribution function, CDF(r)')
ax2.set_ylabel('probability distribution function from cdf, PDF(r)')
# legend_line = ax.legend(handles=leg, loc = 'upper left', numpoints=1)
ax1.axis([t_min, t_max, -0.1, 1])
ax2.axis([t_min, t_max, -0.1, max(diff_cheby)+0.1])
# ax1.set_zorder(ax2.get_zorder() + ax1.
# ax2.set_zorder(ax2.get_zorder())

# following zorder setting is visibility of axis is prefer for the first axis, ax1, rather than ax2
ax1.set_zorder(ax2.get_zorder() + 1)
# However, after changing, the canvas make fill with white blank. That make problem of visibility of the second axis, that is ax2. So, the visible for ax1 is disabled
ax1.patch.set_visible(False)

# Note that the legend should set with ax1 rather than ax2 because of the filling white blank.
ax1.legend(handles=leg_origin, loc = 'upper left', numpoints=1, ncol=1, prop=fontP)

plt.savefig(sys.argv[2])
# plt.show()

# executive parallel computing command:
# parallel python plot_tmp.py {2} cdf_figures/{2}_K{1}.pdf {1} 0.5 1.05 ::: $(seq 8 1 14) ::: $(ls *.info)
# the previous comment is taking all the file information as {2} variable and the degree of Chebyshev-polynomial as {1}. Then output automatically recorded on cdf_figures folder with {2}_K{1}.pdf filename.
