from numpy import *
import scipy.linalg as lin

def get_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
    return kd[argmin(abs(kd))] + x

def get_pos(traj, ts, i, N_dimension):
    index = 2*N_dimension*i + 1
    pos = zeros(N_dimension)
    for k in range(N_dimension):
        pos[k] = traj[ts, index + k]
    return pos

def get_rel_vec(traj, ts, i, j, N_dimension, box_dimension):
    p_i = get_pos(traj, ts, i, N_dimension)
    p_j = get_pos(traj, ts, j, N_dimension)
    for k in range(N_dimension):
        p_j[k] = get_minimum_distance_k_from_x(p_i[k], p_j[k], box_dimension)
    return p_j - p_i

def get_ddf(traj, ts, Np, N_dimension, box_dimension, cut_ratio):
    ddf = []
    for t in ts:
        for i in range(Np-1):
            for j in range(i+1, Np):
                d = lin.norm(get_rel_vec(traj, t, i, j, N_dimension, box_dimension))
                if d<cut_ratio*box_dimension:
                    ddf.append(d)
    return ddf

def get_ddf_angle(traj, ts, Np, N_dimension, box_dimension, cut_ratio, vec_u):
    ddf = []
    angles = []
    for t in ts:
        if t%10 == 0:
            print 'processing %d out of %d'%(t, ts[-1])
        for i in range(Np-1):
            for j in range(i+1, Np):
                tmp_vec = get_rel_vec(traj, t, i, j, N_dimension, box_dimension)
                d = lin.norm(tmp_vec)
                if d<cut_ratio*box_dimension:
                    ddf.append(d)
                    angles.append(arccos(dot(tmp_vec/d, vec_u)))
    return ddf, angles

def get_xy_df(traj, ts, Np, N_dimension, box_dimension, cut_ratio):
    xdf = []
    ydf = []
    cut_distance = box_dimension*cut_ratio
    for t in ts:
        if t%10 == 0:
            print 'processing %d out of %d'%(t, ts[-1])
        
        for i in range(Np-1):
            for j in range(i+1, Np):
                # d = lin.norm(tmp_vec)
                tmp_vec = get_rel_vec(traj, t, i, j, N_dimension, box_dimension)
                
                if tmp_vec[0] < cut_distance and tmp_vec[1] < cut_distance:
                    xdf.append(tmp_vec[0])
                    ydf.append(tmp_vec[1])
    return xdf, ydf

def get_rdf(traj, ts, dr, Np, N_dimension, box_dimension, cut_ratio):
    Nr = int(cut_ratio*box_dimension/dr)
    rdf = zeros([Nr, 3])
    rdf[:,0] = arange(0, cut_ratio*box_dimension, dr)
    ddf = get_ddf(traj, ts, Np, N_dimension, box_dimension, cut_ratio)
    for r in ddf:
        rdf[int(r/dr), 1] += 1
    if (N_dimension == 3):
        rdf[:,2] = rdf[:,1]/(4.*pi*rdf[:,0]**2.0)
        rho_local = size(ddf)/(0.5*(Np-1)*(4./3.)*pi*(cut_ratio*box_dimension)**3.0)
    elif (N_dimension == 2):
        rdf[:,2] = rdf[:,1]/(2.*pi*rdf[:,0])
        rho_local = size(ddf)/(0.5*(Np-1)*pi*(0.5*box_dimension)**2.0)
    normal_fac = dr*size(ts)
    rdf[:,2] /= rho_local*normal_fac
    return rdf

# def get_rdf_angle(traj, ts, dr, Np, N_dimension, box_dimension, cut_ratio, vec_u):
#     Nr = int(cut_ratio*box_dimension/dr)
#     rdf = zeros([Nr, 3])
#     rdf[:,0] = arange(0, cut_ratio*box_dimension, dr)
#     ddf, angles = get_ddf_angle(traj, ts, Np, N_dimension, box_dimension, cut_ratio, vec_u)
#     for r in ddf:
#         rdf[int(r/dr), 1] += 1
#     if (N_dimension == 3):
#         rdf[:,2] = rdf[:,1]/(4.*pi*rdf[:,0]**2.0)
#         rho_local = size(ddf)/(0.5*(Np-1)*(4./3.)*pi*(cut_ratio*box_dimension)**3.0)
#     elif (N_dimension == 2):
#         rdf[:,2] = rdf[:,1]/(2.*pi*rdf[:,0])
#         rho_local = size(ddf)/(0.5*(Np-1)*pi*(0.5*box_dimension)**2.0)
#     normal_fac = dr*size(ts)
#     rdf[:,2] /= rho_local*normal_fac
#     return rdf



def get_rdf_ref(traj, ts, dr, Np, N_dimension, box_dimension, cut_ratio):
    Nr = int(cut_ratio*box_dimension/dr)
    rdf = zeros([Nr, 3])
    rdf[:,0] = arange(0, cut_ratio*box_dimension, dr)
    # ddf = get_ddf(traj, ts, Np, N_dimension, box_dimension, cut_ratio)
    N_tot = 0
    Nt = size(ts)
    for t in ts:
        if t%10 == 0:
            print 'processing %d out of %d'%(t, ts[-1])
        
        for i in range(Np-1):
            for j in range(i+1, Np):
                d = lin.norm(get_rel_vec(traj, t, i, j, N_dimension, box_dimension))
                if d<cut_ratio*box_dimension:
                    rdf[int(d/dr), 1] += 1
                    N_tot += 1
    # for r in ddf:
    #     rdf[int(r/dr), 1] += 1
    if (N_dimension == 3):
        Vr = (4./3.)*pi*((rdf[:,0]+dr)**3.0 - rdf[:,0]**3.0)
        Vrmax = (4./3.)*pi*(cut_ratio*box_dimension)**3.0
        rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    elif (N_dimension == 2):
        Vr = pi*((rdf[:,0]+dr)**2.0 - rdf[:,0]**2.0)
        Vrmax = pi*(cut_ratio*box_dimension)**2.0
        rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    rdf[:,2] = rdf[:,1]/(Vr*0.5*(Np-1)*Nt*rho_local)
    rdf[0, 2] = 0. #removing the first term as zero because it is given by nan (division was 0)
    return rdf, rho_local



def get_rdf_from_ddf(ddf, dr, Np, Nt, N_dimension, box_dimension, cut_ratio):
    Nr = int(cut_ratio*box_dimension/dr) + 1
    rdf = zeros([Nr, 3])
    rdf[:,0] = arange(0, cut_ratio*box_dimension + dr, dr)
    # ddf = get_ddf(traj, ts, Np, N_dimension, box_dimension, cut_ratio)
    N_tot = size(ddf)
    # Nt = size(ts)
    for r in ddf:
        if r < box_dimension*cut_ratio:
            rdf[int(r/dr), 1] += 1
    if (N_dimension == 3):
        Vr = (4./3.)*pi*((rdf[:,0]+dr)**3.0 - rdf[:,0]**3.0)
        Vrmax = (4./3.)*pi*(cut_ratio*box_dimension)**3.0
        rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    elif (N_dimension == 2):
        Vr = pi*((rdf[:,0]+dr)**2.0 - rdf[:,0]**2.0)
        Vrmax = pi*(cut_ratio*box_dimension)**2.0
        rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    rdf[:,2] = rdf[:,1]/(Vr*0.5*(Np-1)*Nt*rho_local)
    rdf[0, 2] = 0. #removing the first term as zero because it is given by nan (division was 0)
    return rdf, rho_local


def get_rdf_from_traj_angle(traj, ts, dr, Np, Nt, N_dimension, box_dimension, cut_ratio, vec_u, dtheta):
    # print cut_ratio, box_dimension, dr
    Nr = int(cut_ratio*box_dimension/dr) + 1
    N_theta = int(2*pi/dtheta) + 1
    rdf = zeros([1 + Nr, 1 + N_theta]) # 1 + for Nr means describing theta index. 1 + for N_theta means describing r index
    rdf[1:,0] = arange(0, cut_ratio*box_dimension + dr, dr)
    rdf[0, 1:] = arange(0, 2*pi, dtheta)
    
    ddf, angles = get_ddf_angle(traj, ts, Np, N_dimension, box_dimension, cut_ratio, vec_u)
    N_tot = size(ddf)
    for i in range(size(ddf)):
        r = ddf[i]
        t = angles[i]
        if r < box_dimension*cut_ratio:
            rdf[ 1 + int(r/dr), 1 + int(t/dtheta)] += 1
            rdf[ 1 + int(r/dr), 1 + int((t + pi)/dtheta)] += 1 # making symmetric
    Vrmax = (4./3.)*pi*(cut_ratio*box_dimension)**3.0
    rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    for i in range(1, 1 + Nr):
        for j in range(1, 1 + N_theta):
            V_ij = (2./3.)*pi*((rdf[i, 0] + dr)**3.0 - rdf[i, 0]**3.0)*(cos(rdf[0, j]) - cos(rdf[0, j] + dtheta))
            rdf[i, j] /= V_ij * 0.5*(Np-1)*Nt*rho_local
    return rdf, rho_local

def get_xydf_from_traj(traj, ts, dr, Np, Nt, N_dimension, box_dimension, cut_ratio):
    # print cut_ratio, box_dimension, dr
    Nr = int(cut_ratio*box_dimension/dr)*2 + 1
    # N_theta = int(2*pi/dtheta) + 1
    rdf = zeros([1 + Nr, 1 + Nr]) # 1 + for Nr means describing theta index. 1 + for N_theta means describing r index
    base_coord = cut_ratio*box_dimension
    rdf[1:,0] = arange(-base_coord, base_coord + dr, dr)
    rdf[0, 1:] = arange(-base_coord, base_coord + dr, dr)

    cut_distance = base_coord
    N_tot = 0
    for t in ts:
        if t%10 == 0:
            print 'processing %d out of %d'%(t, ts[-1])
        
        for i in range(Np-1):
            for j in range(i+1, Np):
                # d = lin.norm(tmp_vec)
                tmp_vec = get_rel_vec(traj, t, i, j, N_dimension, box_dimension)
                
                if tmp_vec[0] < cut_distance and tmp_vec[1] < cut_distance:
                    rdf[ 1 + int((base_coord + tmp_vec[0])/dr), 1 + int((base_coord + tmp_vec[1])/dr)] += 1
                    N_tot += 1
                    # xdf.append(tmp_vec[0])
                    # ydf.append(tmp_vec[1])
    
    # xdf, ydf = get_xy_df(traj, ts, Np, N_dimension, box_dimension, cut_ratio)
    # ddf, angles = get_ddf_angle(traj, ts, Np, N_dimension, box_dimension, cut_ratio, vec_u)
    # N_tot = size(xdf)
    print N_tot
    # for x in xdf:
    #     for y in ydf:
    #         rdf[ 1 + int(x/dr), 1 + int(y/dr)] += 1            
    # for i in range(size(xdf)):
    #     x = xdf[i]
    #     for j in range(size(ydf)):
    #         y = ydf[j]

        # rdf[ 1 + int(/dr), 1 + int((t + pi)/dtheta)] += 1 # making symmetric
    # Vrmax = (4./3.)*pi*(cut_ratio*box_dimension)**3.0
    Vrmax = (box_dimension*cut_ratio)**3.0
    rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    for i in range(1, 1 + Nr):
        for j in range(1, 1 + Nr):
            V_xy = box_dimension*dr**2.0
            rdf[i, j] /= V_xy * 0.5 *(Np - 1)*Nt * rho_local
            # V_ij = (2./3.)*pi*((rdf[i, 0] + dr)**3.0 - rdf[i, 0]**3.0)*(cos(rdf[0, j]) - cos(rdf[0, j] + dtheta))
            # rdf[i, j] /= V_ij * 0.5*(Np-1)*Nt*rho_local
    return rdf, rho_local


def get_rdf_from_rpdist(rpdist, dr, Np, Nt, N_dimension, box_dimension, cut_ratio):
    # Nr = int(cut_ratio*box_dimension/dr) + 1
    r_arr = arange(0, cut_ratio*box_dimension + dr, dr)
    Nr = size(r_arr)
    rdf = zeros([Nr, 3])
    rdf[:,0] = r_arr
    
    # ddf = get_ddf(traj, ts, Np, N_dimension, box_dimension, cut_ratio)
    N_tot = size(rpdist)
    # Nt = size(ts)
    cnt = 0
    for r in rpdist:
        if cnt % int(Nr/100) == 0:
            print '%d out of %d'%(cnt, Nr)
        cnt += 1

        if r < box_dimension*cut_ratio:
            rdf[int(r/dr), 1] += 1
    if (N_dimension == 3):
        Vr = (4./3.)*pi*((rdf[:,0]+dr)**3.0 - rdf[:,0]**3.0)
        Vrmax = (4./3.)*pi*(cut_ratio*box_dimension)**3.0
        rho_local = N_tot/(Nt*(Np-1)*Vrmax)
    elif (N_dimension == 2):
        Vr = pi*((rdf[:,0]+dr)**2.0 - rdf[:,0]**2.0)
        Vrmax = pi*(cut_ratio*box_dimension)**2.0
        rho_local = N_tot/(Nt*0.5*(Np-1)*Vrmax)
    # rdf[:,2] = rdf[:,1]/(Vr*0.5*(Np-1)*Nt*rho_local)
    rho_global = Np/(box_dimension)**3.0
    rdf[:,2] = rdf[:,1]/(Vr*0.5*(Np-1)*Nt*rho_global)
    rdf[0, 2] = 0. #removing the first term as zero because it is given by nan (division was 0)
    return rdf, rho_local

