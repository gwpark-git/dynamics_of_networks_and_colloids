
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def ticks(y, pos):
    return r'$e^{:.0f}$'.format(log(y))

def balance_flux(x_arr, z):
    re = zeros(size(x_arr))
    for i,x in enumerate(x_arr):
        re[i] = (x + z*(1.-x))/(x*(z-1.))
    return re
    # return (x + z(1.-x))/(x*(z-1.))
    # return (z*x + 1.)/((z - 1.)*(1.-x))

x = linspace(0, 1, 100)
Np = 640
Z = Np
y = balance_flux(x, Z)

ref_unity = asarray([[1e-2, 1],
                     [1, 1]])
ref_half = asarray([[0.5, exp(-7)],
                    [0.5, exp(5)]])
# plt.cla()
# plt.clf()
plt.ion()
plt.clf()
# ax.cla()
# ax.clf()
fig, ax = plt.subplots()
ax.loglog(x, y, 'b-', label = 'roof_flux/bridge_flux', basex=e, basey=e)
ax.loglog(ref_unity[:,0], ref_unity[:,1], 'k--', linewidth=2, label = 'ref_unity', basex=e, basey=e)
ax.loglog(ref_half[:,0], ref_half[:,1], 'r--', linewidth=2, label = 'exp(-%4.3f) = 0.5'%(-log(0.5)), basex=e, basey=e)

# ax.set_xscale('ln')
# ax.set_yscale('ln')
# ax.semilogy(x, y, 'b-', label = 'roof_flux/bridge_flux', basey=e)
# ax.semilogy(ref_unity[:,0], ref_unity[:,1], 'k--', linewidth=2, label = 'ref_unity', basey=e)

# ax.plot(x, y, 'b-', label='roof_flux/bridge_flux')
ax.fill_between(x, y, 1, where = y > 1, facecolor = 'blue', alpha = 0.2, interpolate = True)
ax.fill_between(x, y, 1, where = y < 1, facecolor = 'red' , alpha = 0.2, interpolate = True)

ax.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))

# ax.fill_between(log(x), log(y), log(1), where= log(y) > log(1), facecolor = 'blue', interpolate = True, alpha=0.2)
# plt.fill_between(x, y, 1, where= y < 1 , facecolor = 'red', alpha=0.2)
plt.xlabel('probability for energy barrier: exp(-Eb)')
plt.ylabel('ratio between flux for roof and association')
# ax.axis([log(x[0]), log(x[-2]), log(y[0]), log(y[-2])])

# ax.loglog(exp(10), exp(10), 'bs', markersize=10, markerfacecolor='blue', markeredgecolor='none', alpha=0.2, label='prefer for roof', basex=e, basey=e)
# ax.loglog(exp(10), exp(10), 'rs', markersize=10, markerfacecolor='red' , markeredgecolor='none', alpha=0.2, label='prefer for bridge', basex=e, basey=e)

ax.annotate('flux_roof > flux_bridge', xy = (exp(-4), exp(1.9)), size='large')
ax.annotate('flux_roof < flux_bridge', xy = (exp(-0.1), exp(-1)), xytext=(exp(-2), exp(-2)), arrowprops=dict(arrowstyle='->'), size='large')
plt.legend(loc = 'lower left', numpoints=1)
plt.axis([exp(-4.5), exp(0), y[-2], y[1]])
plt.grid()
plt.show()
