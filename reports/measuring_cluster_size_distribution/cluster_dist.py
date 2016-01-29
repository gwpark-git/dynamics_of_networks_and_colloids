
from numpy import *

def check_ident(pos, index, target, Ld):
    Nd = shape(pos)[1]
    if int(target)==-1:
        return -1
    for k in range(Nd):
        if(ident_minimum_distance_k_from_x(pos[index, k], pos[target, k], Ld) != 0):
            return 1
    return 0

def check_travel_beyond_box(pos, index, target, Ld):
    Nd = shape(pos)[1]
    for k in range(Nd):
        if (ident_minimum_distance_k_from_x(pos[index, k], pos[target, k], Ld) != 0):
            return 1
    return 0


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
            # print 'CK: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, N(queue)=%ld'%(cnt, index, target, order_count, size(queue))
            # print queue

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

def DFS_percolation_edge_restricted_box_test(hash, pos, Ld, stack_component, index=0, order_count=1, cnt=0, IDPC=[], IDPI=[], queue=[], queue_order= []):
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
            # print 'CK: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, N(queue)=%ld'%(cnt, index, target, order_count, size(queue))
            # print queue
            tmp_ID = check_ident(pos, index, target, Ld)
            # if int(target) is -1 or tmp_ID is 1:
            if target == -1:
                if size(queue) == 1:
                    return 0
                else:
                    queue = queue[:-1]
                    queue_order = queue_order[:-1]
                    DFS_percolation_edge_restricted_box_test(hash, pos, Ld, stack_component, queue[-1], queue_order[-1] + 1, cnt, IDPC, IDPI, queue, queue_order)
            elif (target in stack_component) or (tmp_ID == 1):
                if tmp_ID == 1:
                    for id in range(shape(pos[index, :])[0]):
                        ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                        if (int(ident_IDP) is not 0) and ([index, target] not in IDPI):
                            IDPC.append([id, ident_IDP])
                            IDPI.append([index, target])
                DFS_percolation_edge_restricted_box_test(hash, pos, Ld, stack_component, index, order_count + 1, cnt, IDPC, IDPI, queue, queue_order)
            else:
                # for id in range(shape(pos[index, :])[0]):
                #     ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                #     if (int(ident_IDP) is not 0) and ([index, target] not in IDPI) and ([target, index] not in IDPI):
                #         IDPC.append([id, ident_IDP])
                #         IDPI.append([index, target])
                    # print IDP
                # IDP.append(ident_percolation(pos[:, index], pos[:, target]))
                stack_component.append(target)
                queue.append(int(target))
                queue_order.append(order_count)
                DFS_percolation_edge_restricted_box_test(hash, pos, Ld, stack_component, target, const_new_order_count, cnt, IDPC, IDPI, queue, queue_order)
        else:
            if size(queue) == 1:
                return 0
            else:
                queue = queue[:-1]
                queue_order = queue_order[:-1]
                DFS_percolation_edge_restricted_box_test(hash, pos, Ld, stack_component, queue[-1], queue_order[-1] +1, cnt, IDPC, IDPI, queue, queue_order)
    except:
        print 'ERR: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld'%(cnt, index, target, order_count)
    return size(queue)


def ident_over(hash, index, order_count):
    N_cols = shape(hash)[1]
    if order_count >= N_cols:
        return 1
    if int(hash[index, order_count]) is -1:
        return 1
    return 0

def DF_percolation_edge_restricted_box_iter(hash, pos, Ld, stack_component, index=0, order_count=1, cnt=0, IDPC=[], IDPI=[], queue=[], queue_order=[]):
    # cnt += 1
    # const_new_order_count = 1
    # N_cols = shape(hash)[1]
    # if size9stack_component)==0:
    cnt = 0
    const_new_order_count = 1
    N_cols = shape(hash)[1]
    queue.append(int(index))
    queue_order.append(order_count)
    try:
        while(size(queue) > 0):
            if order_count < N_cols:
                target = hash[index, order_count]
                tmp_ID = check_ident(pos, index, target, Ld)
                if target == -1:
                    ident_over = True
                    while(ident_over): # this while phrase is design to capture last queue object also over order_count.                        
                        if size(queue) == 1:
                            return 0
                        else:
                            queue = queue[:-1]; queue_order = queue_order[:-1]
                            index = queue[-1]; order_count = queue_order[-1] + 1;
                            if order_count < N_cols:
                                ident_over = False
                elif (target in stack_component) or (tmp_ID == 1):
                    if tmp_ID == 1:
                        for id in range(shape(pos[index, :])[0]):
                            ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                            if (int(ident_IDP) is not 0) and ([index, target] not in IDPI):
                                IDPC.append([id, ident_IDP])
                                IDPI.append([index, target])
                    index = index; order_count = order_count + 1;
                else:
                    stack_component.append(target)
                    queue.append(int(target)); queue_order.append(order_count)
                    index = target; order_count = const_new_order_count;
            else:
                ident_over = True
                while(ident_over):
                    if size(queue) == 1:
                        return 0
                    else:
                        queue = queue[:-1]; queue_order = queue_order[:-1]
                        index = queue[-1]; order_count = queue_order[-1] + 1;
                        if order_count < N_cols:
                            ident_over = False
    except:
        print 'ERR: CNT=%ld, INDEX=%ld, TARGET=%ld, ORDER_CNT=%ld, size_queue=%ld'%(cnt, index, target, order_count, size(queue))
    return size(queue)

def cluster_edge_DFS_travel_restricted_box_iter(hash, pos, Ld, record_component, index=0, order_count=1, cnt=0, IDPC=[], IDPI=[], stack=[], stack_order=[]):
    cnt = 0; const_new_order_count = 1 # initialisation variables
    N_cols = shape(hash)[1] # limitation for the hash tables
    stack.append(int(index)); stack_order.append(order_count) # initial stacking
    while(size(stack) > 0): # will false when size(stack) is 0 if it is not initial step
        cnt += 1 # temporal counting 
        ident_over_cols = ident_over(hash, index, order_count)
        if ident_over_cols: # in the case that the hash[index, order_count] reaching end (-1 or order_count is over)
            stack = stack[:-1]; stack_order = stack_order[:-1]
            if (size(stack) > 0):
                index = stack[-1]; order_count = stack_order[-1] + 1
        else: # in the case that the hash[index, order_count] is properly defined
            target = hash[index, order_count]
            travel_beyond_box = check_travel_beyond_box(pos, index, target, Ld)
            if (target in record_component) or travel_beyond_box: # when target is in stack stack or travel beyond box boundary
                if travel_beyond_box: 
                    for id in range(shape(pos[index, :])[0]):
                        ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                        if (int(ident_IDP) is not 0) and ([index, target] not in IDPI):
                            IDPC.append([id, ident_IDP])
                            IDPI.append([index, target])
                # when particle is duplicated or travel_beyond_box
                index = index; order_count = order_count + 1;

                # this means it inherit the exist index for bead but increase order_count
                # note that the target for next step is given by hash[index, order_count]
            else: # when the target will stack
                record_component.append(int(target))
                stack.append(int(target)); stack_order.append(order_count) # record element and its order for stack
                index = target; order_count = const_new_order_count; # depth first search
    return size(stack)

def cluster_edge_DFS_travel_direct_image_iter(hash, pos, Ld, record_component, index=0, order_count=1, cnt=0, IDPC=[], IDPI=[], stack=[], stack_order=[]):
    cnt = 0; const_new_order_count = 1 # initialisation variables
    N_cols = shape(hash)[1] # limitation for the hash tables
    stack.append(int(index)); stack_order.append(order_count) # initial stacking
    while(size(stack) > 0): # will false when size(stack) is 0 if it is not initial step
        cnt += 1 # temporal counting 
        ident_over_cols = ident_over(hash, index, order_count)
        if ident_over_cols: # in the case that the hash[index, order_count] reaching end (-1 or order_count is over)
            stack = stack[:-1]; stack_order = stack_order[:-1]
            if (size(stack) > 0):
                index = stack[-1]; order_count = stack_order[-1] + 1
        else: # in the case that the hash[index, order_count] is properly defined
            target = hash[index, order_count]
            travel_beyond_box = check_travel_beyond_box(pos, index, target, Ld)
            if (target in record_component) or travel_beyond_box: # when target is in stack stack or travel beyond box boundary
                if travel_beyond_box: 
                    for id in range(shape(pos[index, :])[0]):
                        ident_IDP = ident_minimum_distance_k_from_x(pos[index, id], pos[target, id], Ld)
                        if (int(ident_IDP) is not 0) and ([index, target] not in IDPI):
                            IDPC.append([id, ident_IDP])
                            IDPI.append([index, target])
                # when particle is duplicated or travel_beyond_box
                index = index; order_count = order_count + 1;

                # this means it inherit the exist index for bead but increase order_count
                # note that the target for next step is given by hash[index, order_count]
            else: # when the target will stack
                record_component.append(int(target))
                stack.append(int(target)); stack_order.append(order_count) # record element and its order for stack
                index = target; order_count = const_new_order_count; # depth first search
    return size(stack)


def print_cluster_info(root_index, IDPC, IDPI, Nd, detail=0):
    # pL = []; pR = []
    print 'analysis including root particle %d' % (root_index)
    flag_left = -1
    flag_right = 1
    
    for n in range(Nd):
        print '\tfor axis %d'%(n)
        p_text_left = '\t  travel left pair: '
        p_text_right = '\t  travel right pair: '
        cnt_left = 0
        cnt_right = 0
        for i in range(shape(IDPC)[0]):
            if int(IDPC[i][0]) is n:
                if int(IDPC[i][1]) is flag_left:
                    if detail:
                        p_text_left += '(%3d, %3d), '%(IDPI[i][0], IDPI[i][1])
                    else:
                        cnt_left += 1
                elif int(IDPC[i][1]) is flag_right:
                    if detail:
                        p_text_right += '(%3d, %3d), '%(IDPI[i][0], IDPI[i][1])
                    else:
                        cnt_right += 1
                else:
                    flag_default = 0
        if not(detail):
            p_text_left += ' %d pairs'%(cnt_left)
            p_text_right += ' %d pairs'%(cnt_right)
        print p_text_left
        print p_text_right
        if cnt_left > 0 and cnt_right > 0:
            print '\tcluster with root %d is percolated along %d axis'%(root_index, n)
    return 0
                    


# sc = []; so = []; IDPC=[]; IDPI=[]
# Np = 40
# Nd = 2
# box_dimension = 10.0
# hash = loadtxt('../../test_step/NP40_LD10P2_C100_T3.hash')
# traj = loadtxt('../../test_step/NP40_LD10P2_C100_T3.traj')
# pos = []
# for k in range(Np):
#     pos.append(traj[k*2*Nd + 1:k*2*Nd + 1 + Nd])
# pos = asarray(pos)

# size = DFS_percolation_edge_restricted_box_test(hash=hash, pos=pos, Ld=box_dimension, stack_component=sc, queue_order=so, index=0, IDPC=IDPC, IDPI=IDPI)
# IDPC=asarray(IDPC)
# IDPI=asarray(IDPI)


# sc = []; so = []; IDP = []
# Nd = 2
# test_hash = loadtxt('test.hash')
# test_pos = loadtxt('test.pos')
# box_dimension = 3.0
# # test_hash = loadtxt('test.hash')
# # test_pos = loadtxt('test.pos')
# size = DFS_percolation_test(hash=test_hash, pos=test_pos, Ld=box_dimension, stack_component=sc, queue_order=so, index=0, IDP=IDP)
