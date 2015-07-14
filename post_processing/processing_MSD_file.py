
from numpy import *
import sys

dat = loadtxt(sys.argv[1])
ND = 2
NP = 100
Nt = shape(dat)[0]

tr_RR = zeros([Nt, 2])
tr_RR[:,0] = dat[:,0]
for i in range(NP):
    for k in range(ND):
        index_Rik = 2*ND*i + 1 + k
        tr_RR[:,1] += (dat[:, index_Rik] - dat[0, index_Rik])**2.
tr_RR[:,1] /= float(NP)
savetxt(sys.argv[2], tr_RR)
