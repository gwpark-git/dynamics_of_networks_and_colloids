
from numpy import *
import matplotlib.pyplot as plt

# def DFS_single_root(hash, stack_component, stack_order, index=0, order_count=1, bp_index=0, cnt=0):
#     cnt += 1
#     const_new_order_count = 1
#     N_cols = shape(hash)[1]
#     print cnt, index, order_count, bp_index, stack_component, stack_order
#     if size(stack_component)==0:
#         order_count = const_new_order_count
#         bp_index = 0
#         stack_component.append(index)
#         stack_order.append(order_count)
#     if order_count < N_cols:
#         target = hash[index, order_count]
#         if int(target) is -1:
#             if bp_index < 1:
#                 return 0
#             else:
#                 # print stack_component[bp_index-1], stack_order[bp_index-1]
#                 DFS_single_root(hash, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] + 1, bp_index - 1, cnt)
#         elif target in stack_component:
#             DFS_single_root(hash, stack_component, stack_order, index, order_count + 1, bp_index, cnt)
#         else:
#             stack_component.append(target)
#             stack_order.append(order_count)
#             # bp_index += 1
#             DFS_single_root(hash, stack_component, stack_order, target, const_new_order_count, bp_index + 1, cnt)
#     else:
#         if bp_index < 1:
#             return 0
#         else:
#             DFS_single_root(hash, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] +1, bp_index - 1, cnt)
        
#     return size(stack_component)

def ident_minimum_distance_k_from_x(x, k, box_dimension):
    kd = asarray([k-box_dimension - x, k-x, k+box_dimension-x])
    return argmin(abs(kd)) - 1 # will return [-1, 0, +1]

# def ident_percolation(x, k, box_dimension):
#     return [ident_minimum_distance_k_from_x

# def ident_percolation(p_idim, p_idim):
#     kd = asarray([

def DFS_percolation_test(hash, pos, Ld, stack_component, stack_order, index=0, order_count=1, bp_index=0, cnt=0, IDP=[]):
    try:
        cnt += 1
        const_new_order_count = 1
        N_cols = shape(hash)[1]
        if size(stack_component)==0:
            order_count = const_new_order_count
            bp_index = 0
            stack_component.append(index)
            stack_order.append(order_count)
        if order_count < N_cols:
            target = hash[index, order_count]
            if int(target) is -1:
                if bp_index < 1:
                    return 0
                else:
                    # print stack_component[bp_index-1], stack_order[bp_index-1]
                    DFS_percolation_test(hash, pos, Ld, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] + 1, bp_index - 1, cnt, IDP)
            elif target in stack_component:
                DFS_percolation_test(hash, pos, Ld, stack_component, stack_order, index, order_count + 1, bp_index, cnt, IDP)
            else:
                for id in range(shape(pos[index, :])[0]):
                    ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                    if int(ident_IDP) is not 0:
                        IDP.append([id, ident_IDP, index, target])
                    # print IDP
                # IDP.append(ident_percolation(pos[:, index], pos[:, target]))
                stack_component.append(target)
                stack_order.append(order_count)
                # bp_index += 1
                DFS_percolation_test(hash, pos, Ld, stack_component, stack_order, target, const_new_order_count, bp_index + 1, cnt, IDP)
        else:
            if bp_index < 1:
                return 0
            else:
                DFS_percolation_test(hash, pos, Ld, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] +1, bp_index - 1, cnt, IDP)
    except:
        print 'ERR: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, BP=%ld'%(cnt, index, target, order_count, bp_index)
    return size(stack_component)

# def DFS_percolation_edge_test(hash, pos, Ld, stack_component, stack_order, index=0, order_count=1, bp_index=0, cnt=0, IDPC=[], IDPI=[]):
#     try:
#         cnt += 1
#         if cnt > 50:
#             return 0
#         const_new_order_count = 1
#         N_cols = shape(hash)[1]
#         if size(stack_component)==0:
#             order_count = const_new_order_count
#             bp_index = 0
#             stack_component.append(index)
#             stack_order.append(order_count)
#         if order_count < N_cols:
#             target = hash[index, order_count]
#             print 'CK: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, BP=%ld, N(stack)=%ld'%(cnt, index, target, order_count, bp_index, size(stack_component))

#             if int(target) is -1:
#                 if bp_index < 1:
#                     return 0
#                 else:
#                     # print stack_component[bp_index-1], stack_order[bp_index-1]
#                     DFS_percolation_edge_test(hash, pos, Ld, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] + 1, bp_index - 1, cnt, IDPC, IDPI)
#             elif target in stack_component:
#                 for id in range(shape(pos[index, :])[0]):
#                     ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
#                     if (int(ident_IDP) is not 0) and ([index, target] not in IDPI) and ([target, index] not in IDPI):
#                         IDPC.append([id, ident_IDP])
#                         IDPI.append([index, target])
                
#                 DFS_percolation_edge_test(hash, pos, Ld, stack_component, stack_order, index, order_count + 1, bp_index, cnt, IDPC, IDPI)
#             else:
#                 for id in range(shape(pos[index, :])[0]):
#                     ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
#                     if (int(ident_IDP) is not 0) and ([index, target] not in IDPI) and ([target, index] not in IDPI):
#                         IDPC.append([id, ident_IDP])
#                         IDPI.append([index, target])
#                     # print IDP
#                 # IDP.append(ident_percolation(pos[:, index], pos[:, target]))
#                 stack_component.append(target)
#                 stack_order.append(order_count)
#                 # bp_index += 1
#                 DFS_percolation_edge_test(hash, pos, Ld, stack_component, stack_order, target, const_new_order_count, bp_index + 1, cnt, IDPC, IDPI)
#         else:
#             if bp_index < 1:
#                 return 0
#             else:
#                 DFS_percolation_edge_test(hash, pos, Ld, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] +1, bp_index - 1, cnt, IDPC, IDPI)
#     except:
#         print 'ERR: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, BP=%ld'%(cnt, index, target, order_count, bp_index)
#     return size(stack_component)

def DFS_percolation_edge_test(hash, pos, Ld, stack_component, index=0, order_count=1, cnt=0, IDPC=[], IDPI=[], queue=[], queue_order= []):
    try:
        cnt += 1
        # if cnt > 30:
        #     return 0
        const_new_order_count = 1
        N_cols = shape(hash)[1]
        if size(stack_component)==0:
            order_count = const_new_order_count
            stack_component.append(index)
            queue.append(int(index))
            queue_order.append(order_count)
        if order_count < N_cols:
            target = hash[index, order_count]
            print 'CK: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, N(queue)=%ld'%(cnt, index, target, order_count, size(queue))
            print queue

            if int(target) is -1:
                if size(queue) == 1:
                    return 0
                else:
                    queue = queue[:-1]
                    queue_order = queue_order[:-1]
                    DFS_percolation_edge_test(hash, pos, Ld, stack_component, queue[-1], queue_order[-1] + 1, cnt, IDPC, IDPI, queue, queue_order)
            elif target in stack_component:
                for id in range(shape(pos[index, :])[0]):
                    ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                    if (int(ident_IDP) is not 0) and ([index, target] not in IDPI) and ([target, index] not in IDPI):
                        IDPC.append([id, ident_IDP])
                        IDPI.append([index, target])
                
                DFS_percolation_edge_test(hash, pos, Ld, stack_component, index, order_count + 1, cnt, IDPC, IDPI, queue, queue_order)
            else:
                for id in range(shape(pos[index, :])[0]):
                    ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                    if (int(ident_IDP) is not 0) and ([index, target] not in IDPI) and ([target, index] not in IDPI):
                        IDPC.append([id, ident_IDP])
                        IDPI.append([index, target])
                    # print IDP
                # IDP.append(ident_percolation(pos[:, index], pos[:, target]))
                stack_component.append(target)
                queue.append(int(target))
                queue_order.append(order_count)
                DFS_percolation_edge_test(hash, pos, Ld, stack_component, target, const_new_order_count, cnt, IDPC, IDPI, queue, queue_order)
        else:
            if size(queue) == 1:
                return 0
            else:
                queue = queue[:-1]
                queue_order = queue_order[:-1]
                DFS_percolation_edge_test(hash, pos, Ld, stack_component, queue[-1], queue_order[-1] +1, cnt, IDPC, IDPI, queue, queue_order)
    except:
        print 'ERR: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld'%(cnt, index, target, order_count)
    return size(queue)

sc = []; so = []; IDPC=[]; IDPI=[]
Np = 40
Nd = 2
box_dimension = 10.0
hash = loadtxt('../../test_step/NP40_LD10P2_C100_T3.hash')
traj = loadtxt('../../test_step/NP40_LD10P2_C100_T3.traj')
pos = []
for k in range(Np):
    pos.append(traj[k*2*Nd + 1:k*2*Nd + 1 + Nd])
pos = asarray(pos)

size = DFS_percolation_edge_test(hash=hash, pos=pos, Ld=box_dimension, stack_component=sc, queue_order=so, index=0, IDPC=IDPC, IDPI=IDPI)
IDPC=asarray(IDPC)
IDPI=asarray(IDPI)


# sc = []; so = []; IDP = []
# Nd = 2
# test_hash = loadtxt('test.hash')
# test_pos = loadtxt('test.pos')
# box_dimension = 3.0
# # test_hash = loadtxt('test.hash')
# # test_pos = loadtxt('test.pos')
# size = DFS_percolation_test(hash=test_hash, pos=test_pos, Ld=box_dimension, stack_component=sc, queue_order=so, index=0, IDP=IDP)
