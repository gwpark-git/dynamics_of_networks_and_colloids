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

def acf_step(f, i):
    Nt = size(f)
    val = 0.
    cnt = 0
    for j in range(0, Nt/2 - 1):
        val += f[j]*f[j+i]
        cnt += 1
    return val/float(cnt)

def acf_fix(f):
    Nt = size(f)
    re = zeros(Nt/2)
    for i in range(Nt/2):
        re[i] = acf_step(f, i)
    return re

def acf_average(dat):
    s_xx, s_xy, s_xz, s_yy, s_yz, s_zz = dat[:,0], dat[:,1], dat[:,2], dat[:, 4], dat[:,5], dat[:,8]
    return (acf_fix(s_xy) + acf_fix(s_yz) + acf_fix(s_xz))/3.

def acf_average_flatten(dat):
    return (acf_fix(dat[:,3]) + acf_fix(dat[:,4]) + acf_fix(dat[:,5]))/3.
    
def acf_section(f, Nt_acf, Nt_uncorr):
    Nt = size(f)
    Nt_init = 0
    re = zeros(Nt_acf)
    cnt = 0
    while(Nt_init <= Nt - 2*Nt_acf):
        try:
            for i in range(Nt_acf):
                # index_i = Nt_init + i # using new initial position            
                for j in range(Nt_acf):
                    # index_j = Nt_init + j # using new initial position
                    re[i] += f[Nt_init + j]*f[Nt_init + j + i]
            Nt_init += Nt_uncorr
            cnt += 1
        except:
            print Nt_init, Nt, Nt_acf, Nt - 2*Nt_acf
            break
    if cnt > 0:
        re /= float(cnt)
    else:
        print Nt_init, Nt - Nt_acf
    return re

def acf_average_section(dat, Nt_acf, Nt_uncorr):
    s_xx, s_xy, s_xz, s_yy, s_yz, s_zz = dat[:,0], dat[:,1], dat[:,2], dat[:, 4], dat[:,5], dat[:,8]
    return (acf_section(s_xy, Nt_acf, Nt_uncorr) + acf_section(s_yz, Nt_acf, Nt_uncorr) + acf_section(s_xz, Nt_acf, Nt_uncorr))/3.
    

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
    return re/float(N/2)

def corr_FFT(x, y):
    Fx = fft.rfft(x)
    tempY = []
    y = list(y)
    while len(y) > 0:
        tempY.append(y.pop())
    Fy = fft.rfft(tempY)
    # note that complex conjugate for 1-dimensional real number is itself 
    return fft.irfft(Fx*Fy)/size(x)
    # Nt = size(x)
    # Fx = fft.fft(x, norm='ortho')
    # Fy = fft.fft(y, norm='ortho')
    # F_Corr_xy = conj(Fx)*Fy
    # return abs(fft.fft(F_Corr_xy, norm='ortho'))
    # Fx = fft.fft(x)
    # Fy = fft.fft(y)
    # from scipy.signal import correlate
    # return correlate(asarray(x), asarray(y), mode='same')

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


if __name__=="__main__":
    import sys
    if size(sys.argv) < 5:
        print 'USAGE: main function for computing acf average over section'
        print 'argv[1] == filename for given virial stress data'
        print 'argv[2] == filename for output acf data'
        print 'argv[3] == time conversion factor (from index to time/tau_0)'
        print 'argv[4] == cutting initial part (extracting out non-equilibriated part)'
        print 'argv[5] == (optional) time interval for block average'
        print 'argv[6] == (optional) time range for autocorrelation'
    else:
        fn_virial_stress = sys.argv[1]
        fn_out_acf = sys.argv[2]
        time_factor = float(sys.argv[3])
        init_cut_time = float(sys.argv[4])
        init_cut = int(init_cut_time/time_factor)
        RF = loadtxt(fn_virial_stress)[init_cut:]        
        if(size(sys.argv) < 7):
            acf_RF = acf_average(RF)
            dat = zeros([size(acf_RF), 2])
            dat[:,0] = time_factor*arange(size(acf_RF))
            dat[:,1] = acf_RF
            savetxt(fn_out_acf, dat)
        else:
            tau_acf = float(argv[6])
            Nt_acf = int(tau_acf/time_factor)
            tau_uncorr = float(argv[5])
            Nt_uncorr = int(tau_uncorr/time_factor)
            dat = zeros([Nt_acf, 2])
            dat[:,0] = time_factor*arange(Nt_acf)
            dat[:,1] = acf_average_section(RF, Nt_acf, Nt_uncorr)
            savetxt(fn_out_acf, dat)


            

