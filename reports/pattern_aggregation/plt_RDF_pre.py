
# import sys
# lib_path = os.path.abspath(os.path.join('..','post_processing'))
# sys.path.append(lib_path)
# from lib_rdf import *

from numpy import *
import matplotlib.pyplot as plt
# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# fontP.set_size('x-small')
def get_peak(r, gr, r_min, r_max):
    arg_r_min = get_arg_from_val(r, r_min)
    arg_r_max = get_arg_from_val(r, r_max)
    return r[argmax(gr[arg_r_min:arg_r_max]) + arg_r_min]

def get_arg_from_val(x, x0):
    return argmin(abs(x-x0))

def get_argmin_x0_x1(x, y, x0, x1):
    arg_x0 = get_arg_from_val(x, x0)
    arg_x1 = get_arg_from_val(x, x1)
    return argmin(y[arg_x0:arg_x1]) + arg_x0

def get_n(r, gr, rho, r_init, r_end_0, r_end_1):
    re = 0.
    re = 0.
    dr = r[1] - r[0]
    arg_rmax = get_argmin_x0_x1(r, gr, r_end_0, r_end_1)
    # print dr, arg_rmax, r0, r1
    for i in range(get_arg_from_val(r, r_init), arg_rmax -1):
    # for i in range(get_arg_from_val(r, r0), get_argmin_x0_x1(r, gr, r0, r1) - 1):
        re += 0.5*dr*(gr[i+1] + gr[i])*r[i]#*r[i]**2.0
        # print i, gr[i], re
        # print '%d: g(r)=%e, val=%e, dr = %e'%(i, gr[i], re, dr)
    return 2*pi*rho*re, r[arg_rmax]
    

# def get_n1(r, gr, r0, r1):
#     return get_n(r, gr, 0., r0, r1)
    # re = 0.
    # dr = r[1] - r[0]
    # arg_rmax = get_argmin_x0_x1(r, gr, r0, r1)
    # # print dr, arg_rmax, r0, r1
    # for i in range(arg_rmax -1):
    # # for i in range(get_arg_from_val(r, r0), get_argmin_x0_x1(r, gr, r0, r1) - 1):
    #     re += 0.5*dr*(gr[i+1] + gr[i])*r[i]**2.0
    #     # print i, gr[i], re
    #     # print '%d: g(r)=%e, val=%e, dr = %e'%(i, gr[i], re, dr)
    # return 2*pi*re

# def get_cn(r, gr, r1_st, r1_end, r2_st, r2_end):
#     re_1 = get_n(r, gr, 0., r1_st, r1_end)
#     re_2 = get_n(r, gr, 
    # re_1 = 0.
    # re_2 = 0.
    # dr = r[1] - r[0]
    # arg_rmax_1 = get_argmin_x0_x1(r, gr, r1_st, r1_end)
    # arg_rmax_2 = get_argmin_x0_x1(r, gr, r2_st, r2_end)


# def cal_Nr(r, gr, rmax):
#     re = 0.
#     for i in range(size(r)-1):
#         if r[i] < rmax:
#             re += 0.5*(r[i+1] - r[i])*(gr[i+1] + gr[i])*r[i]#*(r[i])**2.0
#         else:
#             break
#     return 2*pi*re # 2d
        
        
    # re = 0.
    # return get_n1(r, gr, 0., rmax)

# 2pi*rho*\int_{r_0}^{r_1}r dr = pi*rho*(r_1**2.0 - r_0**2.0)

dat = []
rho = 0.4
RDF_NC05 = loadtxt('RDF_G_NC05_p001.dat')
n1_NC05, cut_n1_NC05 = get_n(RDF_NC05[:,0], RDF_NC05[:,2], rho, 0., 1.0, 2.0)
n2_NC05, cut_n2_NC05 = get_n(RDF_NC05[:,0], RDF_NC05[:,2], rho, cut_n1_NC05, 2.0, 3.0)
n3_NC05, cut_n3_NC05 = get_n(RDF_NC05[:,0], RDF_NC05[:,2], rho, cut_n2_NC05, 3.0, 4.0)
ratio_n1_NC05 = n1_NC05/(pi*rho*(cut_n1_NC05**2.0))
ratio_n2_NC05 = n2_NC05/(pi*rho*(cut_n2_NC05**2.0 - cut_n1_NC05**2.0))
ratio_n3_NC05 = n3_NC05/(pi*rho*(cut_n3_NC05**2.0 - cut_n2_NC05**2.0))
dat.append([05, n1_NC05, n2_NC05, n3_NC05, cut_n1_NC05, cut_n2_NC05, cut_n3_NC05, ratio_n1_NC05, ratio_n2_NC05, ratio_n3_NC05])

RDF_NC06 = loadtxt('RDF_G_NC06_p001_pre.dat')
n1_NC06, cut_n1_NC06 = get_n(RDF_NC06[:,0], RDF_NC06[:,2], rho, 0., 1.0, 2.0)
n2_NC06, cut_n2_NC06 = get_n(RDF_NC06[:,0], RDF_NC06[:,2], rho, cut_n1_NC06, 2.0, 3.0)
n3_NC06, cut_n3_NC06 = get_n(RDF_NC06[:,0], RDF_NC06[:,2], rho, cut_n2_NC06, 3.0, 4.0)
ratio_n1_NC06 = n1_NC06/(pi*rho*(cut_n1_NC06**2.0))
ratio_n2_NC06 = n2_NC06/(pi*rho*(cut_n2_NC06**2.0 - cut_n1_NC06**2.0))
ratio_n3_NC06 = n3_NC06/(pi*rho*(cut_n3_NC06**2.0 - cut_n2_NC06**2.0))
dat.append([10, n1_NC06, n2_NC06, n3_NC06, cut_n1_NC06, cut_n2_NC06, cut_n3_NC06, ratio_n1_NC06, ratio_n2_NC06, ratio_n3_NC06])

RDF_NC07 = loadtxt('RDF_G_NC07_p001_pre.dat')
n1_NC07, cut_n1_NC07 = get_n(RDF_NC07[:,0], RDF_NC07[:,2], rho, 0., 1.0, 2.0)
n2_NC07, cut_n2_NC07 = get_n(RDF_NC07[:,0], RDF_NC07[:,2], rho, cut_n1_NC07, 2.0, 3.0)
n3_NC07, cut_n3_NC07 = get_n(RDF_NC07[:,0], RDF_NC07[:,2], rho, cut_n2_NC07, 3.0, 4.0)
ratio_n1_NC07 = n1_NC07/(pi*rho*(cut_n1_NC07**2.0))
ratio_n2_NC07 = n2_NC07/(pi*rho*(cut_n2_NC07**2.0 - cut_n1_NC07**2.0))
ratio_n3_NC07 = n3_NC07/(pi*rho*(cut_n3_NC07**2.0 - cut_n2_NC07**2.0))
dat.append([15, n1_NC07, n2_NC07, n3_NC07, cut_n1_NC07, cut_n2_NC07, cut_n3_NC07, ratio_n1_NC07, ratio_n2_NC07, ratio_n3_NC07])

RDF_NC08 = loadtxt('RDF_G_NC08_p001_pre.dat')
n1_NC08, cut_n1_NC08 = get_n(RDF_NC08[:,0], RDF_NC08[:,2], rho, 0., 1.0, 2.0)
n2_NC08, cut_n2_NC08 = get_n(RDF_NC08[:,0], RDF_NC08[:,2], rho, cut_n1_NC08, 2.0, 3.0)
n3_NC08, cut_n3_NC08 = get_n(RDF_NC08[:,0], RDF_NC08[:,2], rho, cut_n2_NC08, 3.0, 4.0)
ratio_n1_NC08 = n1_NC08/(pi*rho*(cut_n1_NC08**2.0))
ratio_n2_NC08 = n2_NC08/(pi*rho*(cut_n2_NC08**2.0 - cut_n1_NC08**2.0))
ratio_n3_NC08 = n3_NC08/(pi*rho*(cut_n3_NC08**2.0 - cut_n2_NC08**2.0))
dat.append([20, n1_NC08, n2_NC08, n3_NC08, cut_n1_NC08, cut_n2_NC08, cut_n3_NC08, ratio_n1_NC08, ratio_n2_NC08, ratio_n3_NC08])

RDF_NC09 = loadtxt('RDF_G_NC09_p001_pre.dat')
n1_NC09, cut_n1_NC09 = get_n(RDF_NC09[:,0], RDF_NC09[:,2], rho, 0., 1.0, 2.0)
n2_NC09, cut_n2_NC09 = get_n(RDF_NC09[:,0], RDF_NC09[:,2], rho, cut_n1_NC09, 2.0, 3.0)
n3_NC09, cut_n3_NC09 = get_n(RDF_NC09[:,0], RDF_NC09[:,2], rho, cut_n2_NC09, 3.0, 4.0)
ratio_n1_NC09 = n1_NC09/(pi*rho*(cut_n1_NC09**2.0))
ratio_n2_NC09 = n2_NC09/(pi*rho*(cut_n2_NC09**2.0 - cut_n1_NC09**2.0))
ratio_n3_NC09 = n3_NC09/(pi*rho*(cut_n3_NC09**2.0 - cut_n2_NC09**2.0))
dat.append([25, n1_NC09, n2_NC09, n3_NC09, cut_n1_NC09, cut_n2_NC09, cut_n3_NC09, ratio_n1_NC09, ratio_n2_NC09, ratio_n3_NC09])

RDF_NC10 = loadtxt('RDF_G_NC10_p001.dat')
n1_NC10, cut_n1_NC10 = get_n(RDF_NC10[:,0], RDF_NC10[:,2], rho, 0., 1.0, 2.0)
n2_NC10, cut_n2_NC10 = get_n(RDF_NC10[:,0], RDF_NC10[:,2], rho, cut_n1_NC10, 2.0, 3.0)
n3_NC10, cut_n3_NC10 = get_n(RDF_NC10[:,0], RDF_NC10[:,2], rho, cut_n2_NC10, 3.0, 4.0)
ratio_n1_NC10 = n1_NC10/(pi*rho*(cut_n1_NC10**2.0))
ratio_n2_NC10 = n2_NC10/(pi*rho*(cut_n2_NC10**2.0 - cut_n1_NC10**2.0))
ratio_n3_NC10 = n3_NC10/(pi*rho*(cut_n3_NC10**2.0 - cut_n2_NC10**2.0))
dat.append([25, n1_NC10, n2_NC10, n3_NC10, cut_n1_NC10, cut_n2_NC10, cut_n3_NC10, ratio_n1_NC10, ratio_n2_NC10, ratio_n3_NC10])


dat = asarray(dat)
ref_unity = asarray([[0, 1.0],
                     [20, 1.0]])

# ref_cut_n1_X = asarray([cut_n1_NC05, cut_n1_NC10, cut_n1_NC15, cut_n1_NC20, cut_n1_NC25])
# ref_cut_n1_Y = asarray([RDF_NC05[get_arg_from_val(RDF_NC05[:,0], cut_n1_NC05),2],
#                         RDF_NC10[get_arg_from_val(RDF_NC10[:,0], cut_n1_NC10),2],
#                         RDF_NC15[get_arg_from_val(RDF_NC15[:,0], cut_n1_NC15),2],
#                         RDF_NC20[get_arg_from_val(RDF_NC20[:,0], cut_n1_NC20),2],
#                         RDF_NC25[get_arg_from_val(RDF_NC25[:,0], cut_n1_NC25),2]])
# ref_cut_n2_X = asarray([cut_n2_NC05, cut_n2_NC10, cut_n2_NC15, cut_n2_NC20, cut_n2_NC25])
# ref_cut_n2_Y = asarray([RDF_NC05[get_arg_from_val(RDF_NC05[:,0], cut_n2_NC05),2],
#                         RDF_NC10[get_arg_from_val(RDF_NC10[:,0], cut_n2_NC10),2],
#                         RDF_NC15[get_arg_from_val(RDF_NC15[:,0], cut_n2_NC15),2],
#                         RDF_NC20[get_arg_from_val(RDF_NC20[:,0], cut_n2_NC20),2],
#                         RDF_NC25[get_arg_from_val(RDF_NC25[:,0], cut_n2_NC25),2]])
                        

plt.clf()
plt.ion()
# plt.plot(dat[:,0], dat[:,1], 'bo-')
plt.plot(RDF_NC05[:,0], RDF_NC05[:,2], 'b-', linewidth=1, label = 'NC05, R1=%04.2f, R2=%04.2f, R3=%04.2f'%(get_peak(RDF_NC05[:,0], RDF_NC05[:,2], 0, cut_n1_NC05), get_peak(RDF_NC05[:,0], RDF_NC05[:,2], cut_n1_NC05, cut_n2_NC05), get_peak(RDF_NC05[:,0], RDF_NC05[:,2], cut_n2_NC05, cut_n3_NC05)))
plt.plot(cut_n1_NC05, RDF_NC05[get_arg_from_val(RDF_NC05[:,0], cut_n1_NC05),2], 'b|', markersize=15, label='  CN(0, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC05, n1_NC05, cut_n1_NC05, cut_n2_NC05, n2_NC05, cut_n2_NC05, cut_n3_NC05, n3_NC05))
plt.plot(cut_n2_NC05, RDF_NC05[get_arg_from_val(RDF_NC05[:,0], cut_n2_NC05),2], 'b|', markersize=15, label='RCN(0, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC05, ratio_n1_NC05, cut_n1_NC05, cut_n2_NC05, ratio_n2_NC05, cut_n2_NC05, cut_n3_NC05, ratio_n3_NC05))
plt.plot(cut_n3_NC05, RDF_NC05[get_arg_from_val(RDF_NC05[:,0], cut_n3_NC05),2], 'b|', markersize=15)

plt.plot(RDF_NC06[:,0], RDF_NC06[:,2], 'c-', linewidth=1, label = 'NC06, R1=%04.2f, R2=%04.2f, R3=%04.2f'%(get_peak(RDF_NC06[:,0], RDF_NC06[:,2], 0, cut_n1_NC06), get_peak(RDF_NC06[:,0], RDF_NC06[:,2], cut_n1_NC06, cut_n2_NC06), get_peak(RDF_NC06[:,0], RDF_NC06[:,2], cut_n2_NC06, cut_n3_NC06)))
plt.plot(cut_n1_NC06, RDF_NC06[get_arg_from_val(RDF_NC06[:,0], cut_n1_NC06),2], 'c|', markersize=15, label='  CN(0, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC06, n1_NC06, cut_n1_NC06, cut_n2_NC06, n2_NC06, cut_n2_NC06, cut_n3_NC06, n3_NC06))
plt.plot(cut_n2_NC06, RDF_NC06[get_arg_from_val(RDF_NC06[:,0], cut_n2_NC06),2], 'c|', markersize=15, label='RCN(0, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC06, ratio_n1_NC06, cut_n1_NC06, cut_n2_NC06, ratio_n2_NC06, cut_n2_NC06, cut_n3_NC06, ratio_n3_NC06))
plt.plot(cut_n3_NC06, RDF_NC06[get_arg_from_val(RDF_NC06[:,0], cut_n3_NC06),2], 'c|', markersize=15)

plt.plot(RDF_NC07[:,0], RDF_NC07[:,2], 'g-', linewidth=1, label = 'NC07, R1=%04.2f, R2=%04.2f, R3=%04.2f'%(get_peak(RDF_NC07[:,0], RDF_NC07[:,2], 0, cut_n1_NC07), get_peak(RDF_NC07[:,0], RDF_NC07[:,2], cut_n1_NC07, cut_n2_NC07), get_peak(RDF_NC07[:,0], RDF_NC07[:,2], cut_n2_NC07, cut_n3_NC07)))
plt.plot(cut_n1_NC07, RDF_NC07[get_arg_from_val(RDF_NC07[:,0], cut_n1_NC07),2], 'g|', markersize=15, label='  CN(0, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC07, n1_NC07, cut_n1_NC07, cut_n2_NC07, n2_NC07, cut_n2_NC07, cut_n3_NC07, n3_NC07))
plt.plot(cut_n2_NC07, RDF_NC07[get_arg_from_val(RDF_NC07[:,0], cut_n2_NC07),2], 'g|', markersize=15, label='RCN(0, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC07, ratio_n1_NC07, cut_n1_NC07, cut_n2_NC07, ratio_n2_NC07, cut_n2_NC07, cut_n3_NC07, ratio_n3_NC07))
plt.plot(cut_n3_NC07, RDF_NC07[get_arg_from_val(RDF_NC07[:,0], cut_n3_NC07),2], 'g|', markersize=15)


plt.plot(RDF_NC08[:,0], RDF_NC08[:,2], 'y-', linewidth=1, label = 'NC08, R1=%04.2f, R2=%04.2f, R3=%04.2f'%(get_peak(RDF_NC08[:,0], RDF_NC08[:,2], 0, cut_n1_NC08), get_peak(RDF_NC08[:,0], RDF_NC08[:,2], cut_n1_NC08, cut_n2_NC08), get_peak(RDF_NC08[:,0], RDF_NC08[:,2], cut_n2_NC08, cut_n3_NC08)))
plt.plot(cut_n1_NC08, RDF_NC08[get_arg_from_val(RDF_NC08[:,0], cut_n1_NC08),2], 'y|', markersize=15, label='  CN(0, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC08, n1_NC08, cut_n1_NC08, cut_n2_NC08, n2_NC08, cut_n2_NC08, cut_n3_NC08, n3_NC08))
plt.plot(cut_n2_NC08, RDF_NC08[get_arg_from_val(RDF_NC08[:,0], cut_n2_NC08),2], 'y|', markersize=15, label='RCN(0, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC08, ratio_n1_NC08, cut_n1_NC08, cut_n2_NC08, ratio_n2_NC08, cut_n2_NC08, cut_n3_NC08, ratio_n3_NC08))
plt.plot(cut_n3_NC08, RDF_NC08[get_arg_from_val(RDF_NC08[:,0], cut_n3_NC08),2], 'y|', markersize=15)

plt.plot(RDF_NC09[:,0], RDF_NC09[:,2], 'm-', linewidth=1, label = 'NC09, R1=%04.2f, R2=%04.2f, R3=%04.2f'%(get_peak(RDF_NC09[:,0], RDF_NC09[:,2], 0, cut_n1_NC09), get_peak(RDF_NC09[:,0], RDF_NC09[:,2], cut_n1_NC09, cut_n2_NC09), get_peak(RDF_NC09[:,0], RDF_NC09[:,2], cut_n2_NC09, cut_n3_NC09)))
plt.plot(cut_n1_NC09, RDF_NC09[get_arg_from_val(RDF_NC09[:,0], cut_n1_NC09),2], 'm|', markersize=15, label='  CN(0, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC09, n1_NC09, cut_n1_NC09, cut_n2_NC09, n2_NC09, cut_n2_NC09, cut_n3_NC09, n3_NC09))
plt.plot(cut_n2_NC09, RDF_NC09[get_arg_from_val(RDF_NC09[:,0], cut_n2_NC09),2], 'm|', markersize=15, label='RCN(0, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC09, ratio_n1_NC09, cut_n1_NC09, cut_n2_NC09, ratio_n2_NC09, cut_n2_NC09, cut_n3_NC09, ratio_n3_NC09))
plt.plot(cut_n3_NC09, RDF_NC09[get_arg_from_val(RDF_NC09[:,0], cut_n3_NC09),2], 'm|', markersize=15)

plt.plot(RDF_NC10[:,0], RDF_NC10[:,2], 'r-', linewidth=1, label = 'NC10, R1=%04.2f, R2=%04.2f, R3=%04.2f'%(get_peak(RDF_NC10[:,0], RDF_NC10[:,2], 0, cut_n1_NC10), get_peak(RDF_NC10[:,0], RDF_NC10[:,2], cut_n1_NC10, cut_n2_NC10), get_peak(RDF_NC10[:,0], RDF_NC10[:,2], cut_n2_NC10, cut_n3_NC10)))
plt.plot(cut_n1_NC10, RDF_NC10[get_arg_from_val(RDF_NC10[:,0], cut_n1_NC10),2], 'r|', markersize=15, label='  CN(0, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f,   CN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC10, n1_NC10, cut_n1_NC10, cut_n2_NC10, n2_NC10, cut_n2_NC10, cut_n3_NC10, n3_NC10))
plt.plot(cut_n2_NC10, RDF_NC10[get_arg_from_val(RDF_NC10[:,0], cut_n2_NC10),2], 'r|', markersize=15, label='RCN(0, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f, RCN(%04.2f, %04.2f)=%04.2f'%(cut_n1_NC10, ratio_n1_NC10, cut_n1_NC10, cut_n2_NC10, ratio_n2_NC10, cut_n2_NC10, cut_n3_NC10, ratio_n3_NC10))
plt.plot(cut_n3_NC10, RDF_NC10[get_arg_from_val(RDF_NC10[:,0], cut_n3_NC10),2], 'r|', markersize=15)


plt.plot(ref_unity[:,0], ref_unity[:,1], 'k:', linewidth=3, label = 'ref. unity')
plt.legend(loc = 'upper right', numpoints=1, fontsize='xx-small')
plt.axis([0, 5, 0, 9])
plt.xlabel('dimensionless radius')
plt.ylabel('radial distribution function, RDF')
# plt.xticks(range(int(max(RDF_NC25[:,0]))))
plt.grid()
# plt.grid(which='minor', 
plt.show()
