from numpy import *
import matplotlib.pyplot as plt

import sys
sys.path.append('../../../post_processing')

from acf_fcn import *     # autocorrelation functions
from read_file import *   # read file functionality with connectivity information

# parameter setting
base_file_name = 'NP1000_LD10P3_C100'
output_path    = 'tracking_individual_chain'
Np             = 1000
Nc             = 10
box_dimension  = 10.0
N_dimension    = 3
alpha          = 1.5

N_tot          = Np*Nc

# R_chain = zeros([N_dimension, N_tot])
# R_chain = []
cnt = 0
with open (base_file_name + '.traj', 'r') as f_traj:
    with open (base_file_name + '.chain', 'r') as f_chain:
        while(cnt < 100):
            R_chain_step = []
            pos = read_traj_step(f_traj, Np, N_dimension)
            str_chain = f_chain.readline().split('\t')[:-1]
            for CE_HEAD in range(N_tot):
                CE_TAIL = CE_HEAD + N_tot
                index_head = int(str_chain[CE_HEAD])
                index_tail = int(str_chain[CE_TAIL])
                R_chain_step.append(mapped_rel_vec_Rij(pos[index_head,:], pos[index_tail,:], box_dimension))
            R_chain.append(R_chain_step)
            cnt += 1


