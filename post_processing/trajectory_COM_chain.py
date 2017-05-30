from numpy import *
            
def read_traj(f_traj, Np, N_dimension):
    tmp_str = f_traj.readline().split('\t')[:-1]
    t = float(tmp_str[0])
    pos = zeros([Np, N_dimension])
    cnt = 1 # note that tmp_str[0] is time
    for i in range(Np):
        for k in range(N_dimension):
            pos[i, k] = float(tmp_str[cnt])
            cnt += 1
        # cnt += 2 # this will ignore velocity information
        cnt += N_dimension
    return pos, t
# import sys
# sys.path.append('../post_processing/file_IO.py')
# from file_IO import *

def get_index(CE_table, index_chain, ID_head_tail):
    # ID_head_tail==0 : selected chain end is head
    # ID_head_tail==1 : selected chain end is tail
    N_CE = size(CE_table)/2 # count number of chains
    return CE_table[index_chain + ID_head_tail*N_CE]

def get_COM_chain(traj_micelle, CE_table, index_chain, N_dimension):
    flag_head = 0; flag_tail = 1;
    re = (traj_micelle[get_index(CE_table, index_chain, flag_head), :] + traj_micelle[get_index(CE_table, index_chain, flag_tail), :])/2.
    return re

import sys

if size(sys.argv) < 5:
    print 'USAGE for extracting COM chain diffusion:'
    print 'argv[1] : base file name (note that all the information will be stored)'
    print 'argv[2] : output file name'
    print 'argv[3] : number of particles'
    print 'argv[4] : number of chain ends per micelle'
else:
    fn_base = sys.argv[1]
    fn_out = sys.argv[2]
    f_out = file(fn_out, 'a')
    fn_chain = '%s.chain'%(fn_base)
    fn_traj = '%s.traj'%(fn_base)
    Np = int(sys.argv[3])
    N_dimension = 3
    N_CE_per_micelle = int(sys.argv[4])
    NCE = N_CE_per_micelle*Np
    NC = NCE/2
    with open(fn_traj, 'r') as f_traj:
        with open(fn_chain, 'r') as f_chain:
            COM_chain = []
            cnt = 0
            while(1):
                try:
                    print cnt
                    pos, t = read_traj(f_traj, Np, N_dimension)
                    CE_table = f_chain.readline().replace('\n','').split('\t')
                    # COM_chain = zeros([NC+1, N_dimension])
                    COM_chain_tmp = zeros([1, NC*N_dimension*2 + 1])
                    COM_chain_tmp[0, 0] = t
                    for i in range(NC):
                        i_st = i*N_dimension*2 + 1
                        COM_chain_tmp[0, i_st:i_st+N_dimension] = get_COM_chain(pos, CE_table, i, N_dimension)
                    # COM_chain.append(COM_chain_tmp)
                    savetxt(f_out, COM_chain_tmp, fmt='%9.7e')
                    cnt += 1
                except:
                    break
    # COM_chain = asarray(COM_chain)
    # savetxt('%s'%(fn_out), COM_chain, fmt='%7f')


    # with open(fn_traj, 'r') as f_traj:
    #     with open(fn_hash, 'r') as f_hash:
    #         with open(fn_weight, 'r') as f_weight:
    #             with open(fn_chain, 'r') as f_chain:
    #                 cnt = 0
    #                 while(1):
    #                     try:
    #                         pos = read_traj(f_traj, Np, N_dimension)
    #                         connectivity = read_connectivity(f_index, f_weight, Np)

