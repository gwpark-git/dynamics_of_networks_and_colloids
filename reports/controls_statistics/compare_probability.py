
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
def UG(r):
    return r**2.0

def FG(r):
    return 0.8*r

def Boltzmann_map(r):
    return exp(FG(r))

def Boltzmann_factor(r, Eb):
    return exp(-(Eb - FG(r)))

def Boltzmann_factor_cut(r, Eb):
    BF = Boltzmann_factor(r, Eb)
    for i,p in enumerate(BF):
        if p > 1.0:
            BF[i] = 1.0
    return BF
    


r = linspace(0, 2, 80)
Eb = 0.5
# BM = Boltzmann_map(r)
P_old = Boltzmann_map(r)
for i,p in enumerate(P_old):
    tmp = 9./(9.+P_old[0])
    P_old[i] /= 9+P_old[i]
    # P_old[i] += tmp
    # print tmp

BF = Boltzmann_factor(r, Eb)
BFC = Boltzmann_factor_cut(r, Eb)
ref_unity = asarray([[0, 1.0],
                     [5, 1.0]])

arg_p = 0
for i,x in enumerate(BFC):
    if abs(x - 1.0) < 0.0001 and arg_p == 0:
        arg_p = i

plt.clf()
plt.ion()
# plt.plot(r, BM, 'b-', label = 'Normalized probability')
plt.plot(r, P_old, 'b-', linewidth=2, label = r'original approach')
# plt.plot(r, BF, 'g-', linewidth=2, label = r'$\exp(-(\tilde{E}_B-\tilde{F}\tilde{l})$')
plt.plot(r, BFC/(BFC + 9*exp(-Eb)), 'r-', linewidth=2, label = 'flux*Metropolis')
# plt.plot([10], [10], 's', markersize=10, markerfacecolor='blue', markeredgecolor='none', alpha = 0.2, label = 'ACCEPT')
# plt.plot([10], [10], 's', markersize=10, markerfacecolor='red', markeredgecolor='none', alpha = 0.2, label = 'REJECT')

# plt.fill_between(r, BFC, 1, where=BFC < 1, facecolor='red', alpha=0.2, interpolate=True)
# plt.fill_between(r, BFC, 0, where=BFC > 0, facecolor='blue', alpha=0.2, interpolate=True)
# plt.annotate(r'$\exp(-\tilde{E}_B)$', xy=(r[0], BFC[0]), xytext=(r[0]+0.1, BFC[0]-0.1), size='medium', arrowprops=dict(facecolor='black', shrink=0.05))

# plt.annotate(r'$\tilde{E}_B = \tilde{F}\tilde{l}$', xy=(r[arg_p], BFC[arg_p]), xytext=(r[arg_p]+0.1, BFC[arg_p]-0.3), size='medium', arrowprops=dict(facecolor='black', shrink=0.05))

plt.plot(ref_unity[:,0], ref_unity[:,1], 'k--')
plt.axis([0, 1, 0, 1.5])
plt.grid()

plt.xlabel('dimensionless distance')
plt.ylabel('transition probability')
plt.legend(loc = 'upper left', numpoints=1, ncol=2, prop=fontP)
plt.show()
