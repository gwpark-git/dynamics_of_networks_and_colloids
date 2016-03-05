from numpy import *

def acf_gro_step(f, i):
    Nt = size(f)
    val = 0.
    for j in range(Nt - 1 -i):
        val += f[j]*f[j+i]
    return val

def acf_gro_diff_sec(f):
    Nt = size(f)
    re = zeros(Nt/2)
    for i in range(Nt/2):
        re[i] = acf_gro_step(f, i)/float(Nt - i)
    return re

def acf_gro_step_time_lag(f, i, time_lag):
    Nt = size(f)
    val = 0.
    for j in range(Nt - 1 - time_lag):
        val += f[j]*f[j+i]
    return val

def acf_gro(f):
    Nt = size(f)
    re = zeros(Nt/2)
    for i in range(Nt/2):
        re[i] = acf_gro_step_time_lag(f, i, Nt/2)/float(Nt/2)
    return re

def acf_gro_step_uncor(f, i, K):
    Nt = size(f)
    val = 0.
    for j in range(int(Nt/(2*K))-1):
        val += f[j*K]*f[j*K + i]
    return val

def acf_gro_BAV(f, K):
    Nt = size(f)
    re = zeros(Nt/2)
    for i in range(Nt/2):
        re[i] = acf_gro_step_uncor(f, i, K)/float(Nt/K)
    return re

# def acf_gro_BAV(f, K):
#     # K is time index for uncorrelate time
#     Nt = size(f)
#     re = zeros(Nt/2)
#     for i in range(Nt/2):

# def acf_gro_BAV(f, Nb):
#     Nt = size(f)
#     Nt_block = Nt/Nb
#     re = zeros(Nt_block/2)
#     for i in range(Nb):
#         re += acf_gro(f[i*Nt_block:(i+1)*Nt_block])
#     return re/float(Nb)

def acf_gro_BAV_fLAG(f, Nb, fLAG):
    Nt = size(f)
    Nt_block = Nt/Nb
    re = zeros(Nt_block/2)
    for i in range(Nb):
        re += acf_gro(f[i*Nt_block:(i+1)*Nt_block])
    return re/float(Nb)


def corr(x, y):
    N = size(x)
    re = zeros(N/2)
    for i in range(N/2):
        re[i] = dot(x[:N/2], y[i:i+N/2])
    return re

def acf_ref(x):
    N = size(x)
    re = zeros(N/2)
    for i in range(N/2):
        for j in range(N/2):
            re[i] += x[j]*x[j+i]
    return re

def acf(x):
    return corr(x, x)

# def acf_tav(x, N_blocks):
#     Nt = size(x)
#     Nt_block = int(Nt/N_blocks)
#     re = zeros(Nt_block)
#     for i in range(N_blocks):
#         # print i*Nt_block, (i+1)*Nt_block
#         re += acf_wot(x[i*Nt_block:(i+1)*Nt_block])
#     return re


def pacf(x, k): #partial autocorrelation with order of k
    N = size(x) - k
    re = zeros(N/2)
    for i in range(N/2):
        re[i] = dot(x[:N/2], x[k + i:k + i+N/2])
    return re

def acf_np(x):
    re = correlate(x, x, 'full')
    return re[re.size/2:]

def acf_np_tav(x, N_blocks):
    Nt = size(x)
    Nt_block = int(Nt/N_blocks)
    re = zeros(Nt_block)
    for i in range(N_blocks):
        # print i*Nt_block, (i+1)*Nt_block
        re += acf_np(x[i*Nt_block:(i+1)*Nt_block])
    return re

def corr_instant_t(x, y, ref_t):
    N = size(x)
    re = zeros(N/2)
    for i in range(N/2):
        re[i] = x[ref_t]*y[i+ref_t]
    return re

def acf_instant_t(x, ref_t):
    return corr_instant_t(x, x, ref_t)

def acf_wot(x):
    N = size(x)
    re = zeros(N)
    for i in range(N):
        re[i] = x[0]*x[i]
    return re

