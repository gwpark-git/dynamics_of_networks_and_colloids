import sys
from numpy import *
import matplotlib.pyplot as plt
from acf_fcn import *


if size(sys.argv) == 3:
    dat = loadtxt(sys.argv[1])

    N = shape(dat)[0]
    D = shape(dat)[1]
    re = zeros([N, 2])
    re[:,0] = dat[:,0]
    for i in range(1, D):
        re[:,1] += acf_wot(dat[:,i])
    savetxt(sys.argv[2], re)

else:
    print 'USAGE:'
    print 'sys.argv[1] == input data'
    print 'sys.argv[2] == acf without time average'
