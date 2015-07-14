from numpy import *
import matplotlib.pyplot as plt

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



# plt.show()
marker_style = ['bo', 'ro', 'go', 'ko', 'co']
for t in range(0,Nt, 100):
    plt.clf()
    for i in range(Np):
        plt.plot(p[t][i,0], p[t][i,1], marker_style[i%5])

    plt.axis('equal')
    plt.axis([0, box_dimension[0], 0, box_dimension[1]])
    # plt.plot(p[t][:,0], p[t][:,1], 'bo')
    # plt.axis([0, box_dimension[0]*10, 0, box_dimension[1]*10])
    # plt.axis([-5, 5, -5, 5])
    # plt.axis('scaled')
    # plt.xticks(arange(0, box_dimension[0], 0.1))
    # plt.yticks(arange(0, box_dimension[1], 0.1))
    plt.grid('on')
    plt.xlabel('x dimension')
    plt.ylabel('y dimension')
    plt.savefig('figures/t%06d.pdf'%(t))
    # plt.draw()

# plt.plot(traj[t, 
