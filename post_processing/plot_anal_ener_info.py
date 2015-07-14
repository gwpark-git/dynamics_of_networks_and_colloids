
from numpy import *
import matplotlib.pyplot as plt

x = sort(asarray(x))
N_x = size(x)
t = []
for i in linspace(0, N_x/2-1, 1000):
    t.append(x[int(i)])

# N_t = size(t)
# y = arange(0, N_t, 1)/float(N_t)
savetxt('cdf_r.dat', t)

# plt.ion()
# plt.plot(t, y, 'b-')
# plt.grid('on')
# plt.show()

