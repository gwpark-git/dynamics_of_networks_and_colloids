

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
import sys
from chebyshev import *

# x = loadtxt('cdf_r.dat')
# y = arange(0, size(x),1)/float(size(x))
x = sort(loadtxt(sys.argv[1])[:,5])
# for cnt in range(size(x)):
#     if x[cnt] > 0.5:
#         break;
# x = x[cnt:]
    
for cnt in range(size(x)):
    if x[cnt] > 1.05:
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


ref_line = asarray([[1., -0.1],
                    [1., 5]])



c = get_polynomial_from_Cheby(tnew, ynew, int(sys.argv[3]))
t_cheby = linspace(min(tnew), max(tnew), 60)
y_cheby = gen_from_polynomial(t_cheby, c)
diff_cheby = diff_gen_from_polynomial(t_cheby, c)

plt.clf()
fig, ax1 = plt.subplots()
leg = []
ax2 = ax1.twinx()
leg.append(ax1.plot(t, y, 'bo', markerfacecolor='white', markeredgecolor='blue', label='CDF(r), data')[0])
# leg.append(ax1.plot(tnew, ynew, 'r-', label='CDF(r), interpolated')[0])
leg.append(ax1.plot(t_cheby, y_cheby, 'r.', label='CDF(r), Chebyshev')[0])
# leg.append(ax2.plot(ref_line[:,0], ref_line[:,1], 'r--', linewidth=3, label='Range for Repulsion')[0])

# leg.append(ax2.plot(tnew, diff, 'k-', label = 'PDF(r)')[0])
# leg.append(ax2.plot(tnew, 4.*exp(-U_ij(tnew)), 'g-', label = 'U(r)')[0])

leg.append(ax2.plot(t_cheby, diff_cheby, 'r-', linewidth=3, label='PDF(r), Chebyshev')[0])

ax1.grid('on')
ax1.set_xlabel('distance, r')
ax1.set_ylabel('cumulative distribution function, CDF(r)')
ax2.set_ylabel('probability distribution function from cdf, PDF(r)')
# legend_line = ax.legend(handles=leg, loc = 'upper left', numpoints=1)
plt.legend(handles=leg, loc = 'upper left', numpoints=1, ncol=1)
ax1.axis([t_min, t_max, -0.1, 1])
ax2.axis([t_min, t_max, -0.1, max(diff_cheby)])
plt.savefig(sys.argv[2])
# plt.show()

