
tmp_hash = zeros([shape(hash)[0], shape(hash)[1]])
for i in range(shape(hash)[0]):
    for j in range(shape(hash)[1]):
        if int(hash[i, j]) is not -1:
            tmp_hash[i, j] = 1.

# # array([[ 391.,   85.,    0.,    0.,    0.],
# #        [ 415.,   85.,    0.,    0.,    0.],
# #        [ 632.,   90.,    0.,    0.,    0.],
# #        [ 541.,  101.,    0.,    0.,    0.],
# #        [  19.,  170.,    0.,    0.,    0.],
# #        [   0.,  242.,    0.,    0.,    0.],
# #        [ 282.,  244.,    0.,    0.,    0.],
# #        [  81.,  254.,    0.,    0.,    0.],
# #        [  28.,  255.,    0.,    0.,    0.],
# #        [ 458.,  267.,    0.,    0.,    0.]])
# root_index = 9
# rc = []
# tmp_IDPC, tmp_IDPI = cluster_edge_DFS_travel_restricted_box_iter(hash = hash, pos = pos, Ld = box_dimension, record_component=rc, index = root_index)
# print_cluster_info(root_index, tmp_IDPC, tmp_IDPI, Nd)
