
from numpy import *
# import matplotlib.pyplot as plt
import scipy.special as sp
# from matplotlib.font_manager import FontProperties
import sys
# fontP = FontProperties()
# fontP.set_size('x-small')

def sample_signal(x_peak, y_peak, cut_ident):
    Np = size(x_peak)
    re_x = [x_peak[0]]; re_y = [y_peak[0]]
    tmp_array_x = []; tmp_array_y = []
    # for i in range(1, Np-1):
    i=1
    while (i < Np -1):
        dy = y_peak[i] - re_y[-1]
        if abs(dy) > cut_ident:
            tmp_array_x = [x_peak[i]]; tmp_array_y = [y_peak[i]]
            cnt = i
            IDENT_temp = 1
            while(IDENT_temp):
                IDENT_temp = 0
                try:
                    abs_dy_tmp = abs(y_peak[cnt+1] - y_peak[cnt])
                    cnt += 1                
                    tmp_array_x.append(x_peak[cnt]); tmp_array_y.append(y_peak[cnt])
                except:
                    abs_dy_tmp = cut_ident + 1
                if abs_dy_tmp < cut_ident: IDENT_temp = 1

                
            if dy > 0:
                index_tmp_arr = argmax(tmp_array_y)
            else:
                index_tmp_arr = argmin(tmp_array_y)
            print 'x', index_tmp_arr, tmp_array_x
            print 'y', index_tmp_arr, tmp_array_y
            re_x.append(tmp_array_x[index_tmp_arr]); re_y.append(tmp_array_y[index_tmp_arr])
            i = cnt
        else:
            i += 1

    return re_x, re_y

def recursive_convex(x, y, max_cut, max_recursion):
    Nd = size(x)
    re_x = []
    re_y = []
    cnt = 0
    for i in range(1, Nd-1):
        dy1 = y[i] - y[i-1]
        dy2 = y[i+1] - y[i]
        if dy1*dy2 < 0 and x[i]<max_cut:
            re_x.append(x[i])
            re_y.append(y[i])
        else:
            cnt += 1
    if cnt == 0 or max_recursion == 0:
        return re_x, re_y
    return recursive_convex(re_x, re_y, max_cut)

def gr2Sq_2d(q, r, gr, rho):
    Nq = size(q)
    Nr = size(r)
    re = zeros(Nq)
    for i in range(Nq):
        for j in range(Nr-1):
            r1 = r[j]; r2 = r[j+1]
            dr = r2-r1
            f1 = r1*sp.j0(r1*q[i])*(gr[j] - 1.)
            f2 = r2*sp.j0(r2*q[i])*(gr[j+1] - 1.)
            re[i] += 0.5*dr*(f1 + f2)  #trapezoidal rule
    return rho*re + 1.0

def gr2Sq_3d(q, r, gr, rho):
    Nq = size(q)
    Nr = size(r)
    re = zeros(Nq)
    for i in range(Nq):
        for j in range(Nr - 1):
            r1 = r[j]; r2 = r[j+1];
            dr = r2 - r1
            f1 = r1*sin(r1*q[i])*(gr[j] - 1.)
            f2 = r2*sin(r2*q[i])*(gr[j+1] - 1.)
            re[i] += 0.5*dr*(f1 + f2)/q[i]
    # return re + 1.0
    return 4.*pi*rho*re + 1.0
            

if size(sys.argv) < 6:
    print 'USAGE:'
    print 'argv[1] == input rdf file name'
    print 'argv[2] == output Sq file name'
    # print 'argv[3] == cutting identification for sampling'
    print 'argv[3] == minimum of log10(q)'
    print 'argv[4] == maximum of log10(q)'
    print 'argv[5] == number density of particles'
    print 'argv[6] == spatial dimension'
    # print 'argv[6] == minimum of log10(Sq) (plot window)'
    # print 'argv[7] == maximum of log10(Sq) (plot window)'
    # print 'argv[8] == given density (rho)'
    # print 'argv[9] == given max cut for sampling of peak'
    # print 'argv[10] == box dimension'
else:

    fn = sys.argv[1]
    fn_out = sys.argv[2]
    q_log10_min = float(sys.argv[3])
    q_log10_max = float(sys.argv[4])
    q = 2.*pi*10**linspace(q_log10_min, q_log10_max, 1000)
    rdf_r11 = loadtxt(fn)
    r11 = rdf_r11[:,0]
    gr11 = rdf_r11[:,2]
    r11_max = r11[argmax(gr11)]
    dr11 = r11[1] - r11[0]
    rho = float(sys.argv[5])
    Nd = int(sys.argv[6])
    if Nd == 2:
        Sq11 = gr2Sq_2d(q, r11, gr11, rho)
    elif Nd == 3:
        Sq11 = gr2Sq_3d(q, r11, gr11, rho)
    else:
        print 'Wrong Condition for spatial dimension'
    re = zeros([size(q), 2])
    re[:,0] = q*r11_max/(2.*pi)
    re[:,1] = Sq11
    savetxt(fn_out, re)
    
    # ref_unity = asarray([[2.*pi*10**-4, 1], [2.*pi*10**3, 1]])
    # ref_peak = asarray([[1, min(Sq11)], [1, max(Sq11)]])

    # sampled_peak_q11, sampled_peak_Sq11 = recursive_convex(q, Sq11, given_max_cut*(2.*pi)/r11_max, 0)
    # sampled_peak_q11, sampled_peak_Sq11 = sample_signal(sampled_peak_q11, sampled_peak_Sq11, cut_ident)

    # min_x = min(q)
    # max_x = max(q)
    # min_y = min(Sq11)
    # max_y = max(Sq11)

    # ref_box_dimension = asarray([[1./box_dimension, 10**Sq_log10_min], [1./box_dimension, 10**Sq_log10_max]])
    # ref_rdf_distance = asarray([[1./r11[-1], 10**Sq_log10_min], [1./r11[-1], 10**Sq_log10_max]])
    
    # plt.clf()
    # plt.ion()
    # ax = plt.subplot(111)
    # ax.loglog(q*r11_max/(2.*pi), Sq11, 'b-', label = 'computed S(q), %s'%(fn))
    # ax.fill_between(q*r11_max/(2.*pi), Sq11, 1, where=Sq11 > 1, facecolor='blue', alpha=0.2, interpolate=True)
    # ax.fill_between(q*r11_max/(2.*pi), Sq11, 1, where=Sq11 < 1, facecolor='red', alpha = 0.2, interpolate=True)
    # ax.loglog(ref_unity[:,0], ref_unity[:,1], 'k--', linewidth=2, label = 'ref. unity')
    # ax.loglog(ref_box_dimension[:,0], ref_box_dimension[:,1], 'b--', linewidth=2, label = 'rec. box dimension (1./%3.2f)'%(box_dimension), alpha=0.5)
    # ax.annotate('rec. box dimension\n~%3.2e'%(1./box_dimension), xy=(1./box_dimension, 0.9*10**Sq_log10_max), size='xx-small')
    # ax.loglog(ref_rdf_distance[:,0], ref_rdf_distance[:,1], 'r--', linewidth=2, label = 'rec. longest r in g(r) (1./%3.2f)'%(r11[-1]), alpha=0.5)
    # ax.annotate('rec. max r in g(r)\n~%3.2e'%(1./r11[-1]), xy=(1./r11[-1], 0.95*10**Sq_log10_max), size='xx-small')
    # ref_Sq0 = asarray([[10**q_log10_min, Sq11[0]], [10**q_log10_max, Sq11[0]]])
    # ax.loglog(ref_Sq0[:,0], ref_Sq0[:,1], 'g--', label = r'$\tilde{\rho}\tilde{\kappa}_{\tilde{T}} \approx$ %5.4f $\Rightarrow \tilde{\kappa}_{\tilde{T}} \approx $ %5.4f'%(Sq11[0], Sq11[0]/rho))
    # # ax.loglog(ref_peak[:,0], ref_peak[:,1], 'r--', linewidth=2, label = r'ref. peak, $q\approx 2\pi/r_{max}$')

    # ax.loglog(asarray(sampled_peak_q11)*r11_max/(2.*pi), asarray(sampled_peak_Sq11), 'r.', label=r'peaks, ($q*r_{max}/2\pi, Sq$), $r\approx 2\pi/(q*r_{max})$')
    # for i in range(size(sampled_peak_q11)):
    #     given_qi = sampled_peak_q11[i]*r11_max/(2.*pi)
    #     ref_peak_q11 = asarray([[given_qi, min_y],[given_qi, max_y]])
    #     ax.annotate('(%3.2e, %3.2e)\nr~%3.2e'%(given_qi, sampled_peak_Sq11[i], 1./given_qi), xy=(given_qi, sampled_peak_Sq11[i]), size='xx-small')
    # rep_temp = asarray([[1, 1], [1, 1]])
    # ax.grid()
    # # ax.grid(which='minor', axis='y')
    # ax.set_xlabel(r'$q*r_{max}/2\pi$')
    # ax.set_ylabel('isotropic structure factor (2D), S(q)')

    # def my_formatter_fcn(x, p):
    #     return '%.2f' %(x*(10**p))
    # y_ticks = 10**linspace(Sq_log10_min, Sq_log10_max, 6)[1:-1]
    # ax.set_yticks(y_ticks)
    # ax.set_yticklabels(['%4.2e'%(y) for y in y_ticks], rotation=90)
    # # ax.ticklabel_format(axis='y',style='sci', useOffset=False)
    # ax.legend(loc = 'lower right', prop=fontP, numpoints=1)
    # ax.axis([10**q_log10_min, 10**q_log10_max, 10**Sq_log10_min, 10**Sq_log10_max])
    # # aspectratio=1.5
    # # ratio_default=(ax.get_xlim()[1] - ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0])
    # # ax.set_aspect(ratio_default*aspectratio)
    # fig = plt.gcf()
    # fig.set_size_inches(12, 6)
    
    # plt.savefig(fn_out, bbox_inches='tight')
    # # plt.show()

