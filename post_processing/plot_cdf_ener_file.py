

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
import sys

# x = loadtxt('cdf_r.dat')
# y = arange(0, size(x),1)/float(size(x))
U_x = sort(loadtxt(sys.argv[1])[:,6])
for cnt in range(size(U_x)):
    if U_x[cnt] > 1.5:
        break;
U_x = U_x[:cnt]

U_t = []
for i in linspace(0,size(U_x)-1, 5000):
    U_t.append(U_x[int(i)])
N_t = size(U_t)
y = arange(0, N_t, 1)/float(N_t)
U_t_min = min(U_t)
U_t_max = max(U_t)
tck = interpolate.splrep(U_t,y,k=3,s=0)
U_tnew = linspace(U_t_min, U_t_max, 5000)
# xnew = arange(min(x), max(x), 0.0001)
ynew = interpolate.splev(U_tnew, tck, der=0)
diff = interpolate.splev(U_tnew, tck, der=1)


# ref_line = asarray([[1., -0.1],
#                     [1., 5]])


plt.clf()
fig, ax1 = plt.subplots()
leg = []
ax2 = ax1.twinx()
leg.append(ax1.plot(U_t, y, 'bo', markerfacecolor='white', markeredgecolor='blue', label='CDF(U), data')[0])
leg.append(ax1.plot(U_tnew, ynew, 'r-', label='CDF(U), interpolated')[0])
leg.append(ax2.plot(U_tnew, diff, 'k-', label = 'PDF(r)')[0])

# leg.append(ax2.plot(t_cheby, diff_cheby, 'g-', label='PDF(r), Chebyshev')[0])
# leg.append(ax2.plot(ref_line[:,0], ref_line[:,1], 'r--', linewidth=3, label='Range for Repulsion')[0])

ax1.grid('on')
ax1.set_xlabel('distance, r')
ax1.set_ylabel('cumulative distribution function for energy, CDF(r)')
ax2.set_ylabel('probability distribution function from cdf, PDF(r)')
# legend_line = ax.legend(handles=leg, loc = 'upper left', numpoints=1)
plt.legend(handles=leg, loc = 'upper left', numpoints=1, ncol=2)
# ax1.axis([t_min, t_max, -0.1, 1])
# ax2.axis([t_min, t_max, -0.1, 5])
plt.savefig(sys.argv[2])
# plt.show()

# c = get_polynomial_from_Cheby(tnew, ynew, 12)
# t_cheby = linspace(min(tnew), max(tnew), 100)
# y_cheby = gen_from_polynomial(t_cheby, c)
# diff_cheby = diff_gen_from_polynomial(t_cheby, c)
