


from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as lin
import sys
sys.path.append('../post_processing')
from read_file import *


if size(sys.argv) < 4:
    print 'USAGE:'
    print 'argv[1] == hash file name'
    print 'argv[2] == out file name'
    print 'argv[3] == Np'
    # print 'argv[4] == N_dimension'
    # print 'argv[5] == box_dimension'
    # print 'argv[6] == lc'
    # print 'argv[7] == alpha'
    # print 'argv[8] == minimum_time_index (DEFAULT 0)'
    # print 'argv[9] == maximum_time_index (DEFAULT -1)'
else:

    fn_base = sys.argv[1]
    Np = int(sys.argv[3])
    NDAS = []
    with open (fn_base, 'r') as f_index:
        tmp_cnt = 0
        while(1):
            try:
                tmp_cnt +=1
                tmp_NDAS = 0
                for i in range(Np):
                    tmp_str_index = f_index.readline().split('\t')[:-1]
                    for j in range(size(tmp_str_index)):
                        if int(tmp_str_index[j]) == -1:
                            tmp_NDAS += j - 1
                            break
                # NDAS.append(size(tmp_str_index) - 1)
                NDAS.append(tmp_NDAS)
                if tmp_cnt % 100 == 0:
                    print tmp_cnt, tmp_NDAS
            except:
                break
    savetxt(sys.argv[2], NDAS)
