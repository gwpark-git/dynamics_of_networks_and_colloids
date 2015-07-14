
from numpy import *
import matplotlib.pyplot as plt

dat = loadtxt('../data_test_longer/100.info')
path_fig = '../data_test_longer/100.info.figure'
N_rows, N_cols = shape(dat)

# N_ref = 1001
st_index = 0
end_index = -1
plt.clf()
line_color=['b', 'r', 'y', 'c', 'k']
cnt = 0
for i in range(1, N_rows):
    if dat[i,0] <> dat[i-1, 0]:
        end_index = i
        distance_part = sort(dat[st_index:end_index, 5])
        Nr = arange(0, size(distance_part), 1)
        plt.clf()
        plt.plot(distance_part, Nr/float(size(distance_part)), 'b-', color=line_color[cnt%size(line_color)], label='steps=%d'%(dat[i-1, 0]))
        plt.legend(loc = 'upper left')
        plt.axis([0.0, 3, 0, 1])
        plt.grid('on')
        plt.xlabel('reduced distance')
        plt.ylabel('Cumulative Distribution')
        plt.savefig(path_fig+'/tmp_fig_%08d.png'%(cnt), dpi=300)
        plt.close()
        st_index = end_index + 1
        cnt += 1
# plt.show()        


