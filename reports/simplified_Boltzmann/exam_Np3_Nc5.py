

from numpy import *
import matplotlib.pyplot as plt
from scipy.linalg import norm

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
Nc = 5
# ((Np, Nc)) = (Np + Nc - 1, Nc) => 7!/(5!2!) = 21
Ns = 21

# w0 = zeros([Ns, Np])
# w1 = zeros([Ns, Np])
# w2 = zeros([Ns, Np])

w = asarray([[5, 0, 0], #0
             [4, 1, 0], [4, 0, 1], #12
             [3, 2, 0], [3, 1, 1], [3, 0, 2], #345
             [2, 3, 0], [2, 2, 1], [2, 1, 2], [2, 0, 3], #6789
             [1, 4, 0], [1, 3, 1], [1, 2, 2], [1, 1, 3], [1, 0, 4], #10-14
             [0, 5, 0], [0, 4, 1], [0, 3, 2], [0, 2, 3], [0, 1, 4], [0, 0, 5]]) #15-20


tmp_P0 = zeros(Ns)
tmp_P1 = zeros(Ns)
tmp_P2 = zeros(Ns)
for i in range(Ns):
    for j in range(Np):
        tmp_P0[i] += w[i,j]*uij(r[0, j])
        tmp_P1[i] += w[i,j]*uij(r[1, j])
        tmp_P2[i] += w[i,j]*uij(r[2, j])
    tmp_P0[i] = exp(-tmp_P0[i])
    tmp_P1[i] = exp(-tmp_P1[i])
    tmp_P2[i] = exp(-tmp_P2[i])

Z0, Z1, Z2 = sum(tmp_P0), sum(tmp_P1), sum(tmp_P2)
P0, P1, P2 = tmp_P0/Z0,   tmp_P1/Z1,   tmp_P2/Z2


plt.clf()
plt.ion()
plt.plot(P0, 'b.-', label='P0, max=(%d, %4.3f)'%(argmax(P0), P0[argmax(P0)]))
plt.plot(P1, 'r.-', label='P1, max=(%d, %4.3f)'%(argmax(P1), P1[argmax(P1)]))
plt.plot(P2, 'g.-', label='P2, max=(%d, %4.3f)'%(argmax(P2), P2[argmax(P2)]))
plt.grid()
plt.legend(loc = 'upper left')
plt.xlabel('index for microstate')
plt.ylabel('probability for microstate')
plt.axis([0, Ns-1, -0.01, 1.0])
plt.show()
