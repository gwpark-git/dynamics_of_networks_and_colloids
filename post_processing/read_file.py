
from numpy import *
import itertools as iter


# class condition:
#     """ the condition will read the input file of the infrastructure source code """
    

#     def __init__(self, 


def read_chunks(f, chunks=100):
    # read files with certain chunks
    # the size of chunk is the second argumnet and default is set with 100
    block = True
    while block:
        block = [f.readline() for i in range(chunks)]
        block = list(filter(None, block))
        yield block

def get_seq_data_int(fn, skip_w = 1, ncols = -1, skiprows = 0):
    # skip_w means skipping frequency
    # ncols means the number of columns
    # skiprows means skipping header number
    with open (fn, 'r') as f:
        if ncols is -1:
            dat = genfromtxt(iter.islice(f, 0, None, skip_w), dtype=int, skip_header = skiprows)
        else:
            c_arr = range(ncols)
            dat = genfromtxt(iter.islice(f,0, None, skip_w), dtype=int, skip_header = skiprows, usecols=c_arr)
    return dat

def get_seq_data_float(fn, skip_w = 1, ncols = -1, skiprows = 0):
    # skip_w means skipping frequency
    # ncols means the number of columns
    # skiprows means skipping header number
    with open (fn, 'r') as f:
        if ncols is -1:
            dat = genfromtxt(iter.islice(f, 0, None, skip_w), dtype=float, skip_header = skiprows)
        else:
            c_arr = range(ncols)
            dat = genfromtxt(iter.islice(f,0, None, skip_w), dtype=float, skip_header = skiprows, usecols=c_arr)
    return dat


def get_block_seq_data_int(fn, nblocks = 1, skip_w = 1, ncols = -1, skiprows = 0):
    with open (fn, 'r') as f:
        dat = []
        for cnt, piece in enumerate(read_chunks(f, nblocks)):
            if cnt <> 0 and cnt%skip_w == 0:
                dat += piece
    if ncols is -1:
        dat = genfromtxt(dat, dtype=int)
    else:
        dat = genfromtxt(dat, dtype=int, usecols = range(ncols))
    return dat
            
                # dat = genfromtxt(piece, dtype=int)
                

def get_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
    return kd[argmin(abs(kd))]+x

def map_minimum_image_Rj_from_Ri(Ri, Rj, box_dimension):
    N_dimension = size(Rj)
    minimum_image_Rj = zeros(N_dimension)
    for i in range(N_dimension):
        minimum_image_Rj[i] = get_minimum_distance_k_from_x(Ri[i], Rj[i], box_dimension)
    return minimum_image_Rj

def mapped_rel_vec_Rij(Ri, Rj, box_dimension):
    mapped_Rj = map_minimum_image_Rj_from_Ri(Ri, Rj, box_dimension)
    return mapped_Rj - Ri

def RR_over_beads(pos, connectivity, box_dimension):
    Np, N_dimension = shape(pos)
    RR = zeros([N_dimension, N_dimension])
    for i in range(Np):
        tmp_index = connectivity[i,i:].nonzero()[0] + i
        #note that nonzero returns tuple. So, the [0] indice will return ndarray
        for j in tmp_index:
            R_j = map_minimum_image_Rj_from_Ri(pos[i,:], pos[j,:], box_dimension)
            R_ij = rel_vec_Rij(pos[i,:], R_j)
            RR += connectivity[i,j]*outer(R_ij, R_ij)
    RR /= float(Np)
    return RR
            
# def read_traj_full(fn_traj, Np, N_dimension):
#     with open (fn_traj, 'r') as f_traj:
#         tmp_str = f_traj.readline().split('\t')[:-1]
#         pos = zeros([Np, N_dimension])
#         cnt = 1 # note that tmp_str[0] is time
#         for i in range(Np):
#             for k in range(N_dimension):
#                 pos[i, k] = float(tmp_str[cnt])
#                 cnt += 1
#             # cnt += 2 # this will ignore velocity information
#             cnt += N_dimension
#     return pos

# def read_connectivity_full(fn_index, fn_weight, Np):
#     with open (fn_index, 'r') as f_index:
#         with open (fn_weight, 'r') as f_weight:
#             connectivity = zeros([Np, Np])
#             str_index_table = []
#             str_weight_table = []
#             for i in range(Np):
#                 str_index_table.append(f_index.readline().split('\t')[:-1])
#                 str_weight_table.append(f_weight.readline().split('\t')[:-1])
#             N_cols = shape(str_index_table)[1]
#             for i in range(Np):
#                 index_i = int(str_index_table[i][0])
#                 for j in range(N_cols):
#                     index_j = int(str_index_table[i][j])
#                     if index_j == -1:
#                         break
#                     else:
#                         connectivity[index_i, index_j] = int(str_weight_table[i][j])
#     return connectivity

def read_traj_step(f_traj, Np, N_dimension):
    tmp_str = f_traj.readline().split('\t')[:-1]
    pos = zeros([Np, N_dimension])
    cnt = 1 # note that tmp_str[0] is time
    for i in range(Np):
        for k in range(N_dimension):
            pos[i, k] = float(tmp_str[cnt])
            cnt += 1
        # cnt += 2 # this will ignore velocity information
        cnt += N_dimension
    return pos

def read_connectivity_step(f_index, f_weight, Np):
    connectivity = zeros([Np, Np])
    str_index_table = []
    str_weight_table = []
    max_N_cols = 1
    for i in range(Np):
        tmp_str_index = f_index.readline().split('\t')[:-1] # -1 is used because of \n
        tmp_str_weight = f_weight.readline().split('\t')[:-1]
        index_i = int(tmp_str_index[0])
        for j in range(size(tmp_str_index)):
            index_j = int(tmp_str_index[j])
            if index_j == -1:
                break
            else:
                connectivity[index_i, index_j] = int(tmp_str_weight[j])
        # if max_N_cols < size(tmp_str_index):
        #     max_N_cols 
        # str_index_table.append(f_index.readline().split('\t')[:-1])
        # str_weight_table.append(f_weight.readline().split('\t')[:-1])
    # N_cols = shape(str_index_table)[1]
    # for i in range(Np):
    #     index_i = int(str_index_table[i][0])
    #     for j in range(N_cols):
    #         index_j = int(str_index_table[i][j])
    #         if index_j == -1:
    #             break
    #         else:
    #             connectivity[index_i, index_j] = int(str_weight_table[i][j])
    return connectivity

                
# def get_block_seq_data_int(fn, nblocks = 1, skip_w = 1, ncols = -1, skiprows = 0):
#     # unlike get_seq_data_int, get_block_seq_data_int will read the block of data in sequencial manner
#     with open (fn, 'r') as f:
#         cnt = skiprows
#         ident = True
#         while ident:
#             block = [f.readline() for i in range(nblocks)]
    # n_arr = range(skiprows, None, skip_w)
    # with open (fn, 'r') as f:
    #     cnt = 0
    #     it_chain = iter.chain(iter.islice(f, cnt), iter.islice(f, nblocks, None))
    #     try:
    #         while(1):
    #             cnt += skip_w
    #             tmp_it_chain = iter.chain(iter.islice(f, cnt), iter.islice(f, nblocks, None))
    #             it_chain = iter.chain(it_chain, tmp_it_chain)
    #     except:
    #         if ncols is -1:
    #             dat = genfromtxt(it_chain, dtype=int, skip_header = skiprows)
    #         else:
    #             dat = genfromtxt(it_chain, dtype=int, skip_header = skiprows, usecols = range(ncols))
    #     return dat
                
        # seed_lines = list(iter.islice(f, 0, None, skip_w))
        # for x in seed_lines:
        #     it.chain(
        # try:
        #     while(1):
        #         tmp = iter.islice(f, 0, None, skip_w)
        #         for i,it in enumerate(tmp):
        #             lines = iter.chain(tmp, iter.islice(f, 
        # lines = iter.islice(f, 0, None, skip_w)
        
        # if ncols is -1:
        #     dat = genfromtxt(iter.chain(iter.islice(f, 0, None, skip_w), iter.islice(f, nblocks, None)), dtype=int, skip_header = skiprows)
        # else:
        #     c_arr = range(ncols)
        #     dat = genfromtxt(iter.chain(iter.islice(f, 0, None, skip_w), iter.islice(f, nblocks, None)), dtype=int, skip_header = skiprows, usecols=c_arr)

    # return dat

# dat = get_data_int('test_step/one_cluster.hash', 15, 10)
