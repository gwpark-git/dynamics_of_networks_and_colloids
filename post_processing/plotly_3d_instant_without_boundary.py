from numpy import *
import pylab as P

import sys

# dat = loadtxt('tmp.traj')
# hash = loadtxt('tmp.hash')
# Np = 3200 # number of particles
# Nd = 3 # spatial dimension

if size(sys.argv) < 7:
    print 'USAGE: 3D plot via plotly package'
    print 'return file in html format'
    print 'argv[1] == trajectory'
    print 'argv[2] == hash'
    print 'argv[3] == path for output htmls'
    print 'argv[4] == spatial dimension'
    print 'argv[5] == number of particles'
    print 'argv[6] == box dimension'
    # print 'argv[7] == df for given trajectory'
else:    
    dat = loadtxt(sys.argv[1])
    # hash = loadtxt(sys.argv[2])
    o_path = sys.argv[3]
    Nd = long(sys.argv[4])
    Np = long(sys.argv[5])
    box_dimension = float(sys.argv[6])
    # df = float(sys.argv[7])
    hash = []

    with open (sys.argv[2], 'r') as f_hash:
        # hash_tmp = zeros([Np, Np])
        cnt = 0
        for line in f_hash:
            hash_tmp = -1*ones(Np)
            tmp_str = line.replace('\t\n','').split('\t')
            for i,h in enumerate(tmp_str):
                hash_tmp[i] = float(h)
                # print tmp_str[:5]
                # print hash_tmp[:5]
            hash.append(hash_tmp)
    # print hash
    hash = asarray(hash)

    from scipy.linalg import norm
    def distance(x_pair, y_pair, z_pair):
        p1 = asarray([x_pair[0], y_pair[0], z_pair[0]])
        p2 = asarray([x_pair[1], y_pair[1], z_pair[1]])
        return norm(p1 - p2)

    ref_unity = asarray([[0, 1], [5, 1]])

    import matplotlib.pyplot as plt
    import plotly
    import plotly.graph_objs as go

    from pylab import rand

    def get_minimum_distance_k_from_x(x, k, box_dimension):
        kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
        return kd[argmin(abs(kd))] + x;

    def check_boundary_jump(x, k, box_dimension):
        kd = asarray([k-box_dimension-x, k-x, k+box_dimension-x])
        if argmin(abs(kd)) <> 1:
            return True
        return False
    
    pos = zeros([Np, Nd])


    color_map = zeros([Np, 3])
    for i in range(Np):
        color_map[i, :] = rand(3)*255
    plotly_color_map = []
    plotly_text = []
    for i in range(Np):
        plotly_color_map.append('rgb(%f, %f, %f)'%(color_map[i,0], color_map[i,1], color_map[i,2]))
        plotly_text.append('index = %d'%(i))


    Nt = shape(dat)[0]

    t_cnt = 0
    # t_arr = [0, 1000, 2000]
    t_arr = range(Nt)
    t_dist = []
    # for t in arange(0, Nt, 10):
    #     t_cnt += 1
    # for t in range(0, 1, 10):
    # for t in t_arr:
    for i in range(Np):
        for k in range(Nd):
            pos[i, k] = dat[2*Nd*i + 1 + k]

    N_dimension = Nd


    tr_particle = go.Scatter3d(
        x=pos[:,0],
        y=pos[:,1],
        z=pos[:,2],
        text=plotly_text,
        legendgroup='particle',
        mode='markers',
        marker=dict(
            size=3,
            color=plotly_color_map,
            line=dict(
                # color='rgba(217, 217, 217, 0.14)',
                # color=color_map,
                color=plotly_color_map,
                width=0.5
                ),
            opacity=0.8
            )
        )

    cnt = 0
    N_cols = shape(hash)[1]
    # tr_association = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=go.Line(color='red'), opacity=0.5)
    trace = [tr_particle]
    hash_st = 0
    d_dist = []
    print hash
    for i in range(Np/2 + (Np%2)): # excluding duplication
        for j in range(1, N_cols): # excluding itself
            # if i > hash[hash_st + i, j]
            # if i <> hash[hash_st + i,j]:
            if hash[hash_st + i,j] != -1:
                # The following scheme is used 'None' in order to generate disconnected line plot
                # Unlike 2d line plot, however, the functionality is not properly working which is weird for me
                # Loosing disconnected line for one trace means we need big overhead to show the plot in web browser.
                # tr_association['x'] += [pos[i, 0], get_minimum_distance_k_from_x(pos[i, 0], pos[hash[i,j], 0], box_dimension), None]
                # tr_association['y'] += [pos[i, 1], get_minimum_distance_k_from_x(pos[i, 1], pos[hash[i,j], 1], box_dimension), None]
                # tr_association['z'] += [pos[i, 2], get_minimum_distance_k_from_x(pos[i, 2], pos[hash[i,j], 2], box_dimension), None]
                if (check_boundary_jump(pos[i, 0], pos[hash[hash_st + i,j], 0], box_dimension) == False and check_boundary_jump(pos[i, 1], pos[hash[hash_st + i,j], 1], box_dimension) == False and check_boundary_jump(pos[i, 2], pos[hash[hash_st + i,j], 2], box_dimension) == False ):
                    x_pair = [pos[i, 0], get_minimum_distance_k_from_x(pos[i, 0], pos[hash[hash_st + i,j], 0], box_dimension)]
                    y_pair = [pos[i, 1], get_minimum_distance_k_from_x(pos[i, 1], pos[hash[hash_st + i,j], 1], box_dimension)]
                    z_pair = [pos[i, 2], get_minimum_distance_k_from_x(pos[i, 2], pos[hash[hash_st + i,j], 2], box_dimension)]
                    d_pair = distance(x_pair, y_pair, z_pair)
                    # d_dist.append([cnt, d_pair])
                    d_dist.append(d_pair)
                    cnt += 1
                    trace.append(go.Scatter3d(x=x_pair, y=y_pair, z=z_pair, legendgroup='bridge', mode='lines', text='p(%d, %d)\nd=%4.3f\n'%(i,j,d_pair)))
            else:
                break

    axis=dict(showbackground=False,
              zeroline=False,
              showgrid=True,
              showticklabels=True,
              title='')


    layout = go.Layout(title=sys.argv[1],
                       margin=go.Margin(t=100),
                       showlegend=False,
                       scene=go.Scene(xaxis=go.XAxis(axis), yaxis=go.YAxis(axis), zaxis=go.ZAxis(axis)))

    fig = go.Figure(data=trace, layout=layout)
    plot_url = plotly.offline.plot(fig, filename='%s/%s.html'%(o_path, sys.argv[1]), auto_open=False)

    d_dist = asarray(d_dist)
    t_dist.append(d_dist)


