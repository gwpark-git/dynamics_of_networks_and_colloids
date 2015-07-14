
from numpy import *
import matplotlib.pyplot as plt
given_dat = loadtxt('../data_test_longer/100.info')
# path_fig = '../data_test_longer/100.info.figure'
N_rows, N_cols = shape(given_dat)

# N_ref = 1001
st_index = 0
end_index = -1
plt.clf()
line_color=['b', 'r', 'y', 'c', 'k']
cnt = 0
for i in range(1, N_rows):
    if given_dat[i,0] <> given_dat[0, 0]:
        cnt = i
        break

x = []
U_x = []
for i in range(cnt, N_rows):
    x.append(given_dat[i, 5])
    U_x.append(given_dat[i, 6])
    
# dat = []
# dat.append(x)
# dat.append(U_x)
# dat = transpose(asarray(dat))
