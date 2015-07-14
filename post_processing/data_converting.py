from numpy import *


traj = loadtxt('tmp.dat')
Nt = shape(traj)[0]
dimension = 2
Np = int((shape(traj)[1] - 1)/(2*dimension)) 
box_dimension = [1.0, 1.0]

p = zeros([Nt, Np, dimension])
v = zeros([Nt, Np, dimension])


for t in range(Nt):
    for i in range(Np):
        for k in range(dimension):
            # print i*dimension*2 + 1 + k, i*dimension*2 + 1 + dimension + k
            p[t, i, k] = traj[t, i*dimension*2 + 1 + k]%box_dimension[k]
            v[t, i, k] = traj[t, i*dimension*2 + 1 + dimension + k]%box_dimension[k]

savetxt(p, 'tmp_position.dat')
savetxt(v, 'tmp_velocity.dat')

