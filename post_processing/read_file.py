
from numpy import *
import itertools as iter

def read_chunks(f, chunks=100):
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
