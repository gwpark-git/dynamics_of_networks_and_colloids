

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate

# x = loadtxt('cdf_r.dat')
# y = arange(0, size(x),1)/float(size(x))
x = sort(loadtxt('tmp_r.dat'))
for cnt in range(size(x)):
    if x[cnt] > 1.5:
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



# c = get_polynomial_from_Cheby(tnew, ynew, 12)
# t_cheby = linspace(min(tnew), max(tnew), 100)
# y_cheby = gen_from_polynomial(t_cheby, c)
# diff_cheby = diff_gen_from_polynomial(t_cheby, c)
