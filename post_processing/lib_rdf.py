
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

def get_rdf(traj, ts, dr, Np, N_dimension, box_dimension, cut_ratio):
    Nr = int(cut_ratio*box_dimension/dr)
    rdf = zeros([Nr, 3])
    rdf[:,0] = arange(0, cut_ratio*box_dimension, dr)
    ddf = get_ddf(traj, ts, Np, N_dimension, box_dimension, cut_ratio)
    for r in ddf:
        rdf[int(r/dr), 1] += 1
    rdf[:,1] /= dr
    if (N_dimension == 3):
        rdf[:,2] = rdf[:,1]/(4.*pi*rdf[:,0]**2.0)
        rdf_local = size(ddf)/(0.5*(Np-1)*(4./3.)*pi*(cut_ratio*box_dimension)**3.0)
    elif (N_dimension == 2):
        rdf[:,2] = rdf[:,1]/(2.*pi*rdf[:,0])
        rho_local = size(ddf)/(Np*pi*(cut_ratio*box_dimension)**2.0)
    print 'rho_local = %f'%(rho_local)
    rdf[:,2] /= rho_local
    return rdf
