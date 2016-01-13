
from numpy import *
import matplotlib.pyplot as plt

def DFS_single_root(hash, stack_component, stack_order, index=0, order_count=1, bp_index=0, cnt=0):
    cnt += 1
    const_new_order_count = 1
    N_cols = shape(hash)[1]
    print cnt, index, order_count, bp_index, stack_component, stack_order
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
                recursive_seq(hash, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] + 1, bp_index - 1, cnt)
        elif target in stack_component:
            recursive_seq(hash, stack_component, stack_order, index, order_count + 1, bp_index, cnt)
        else:
            stack_component.append(target)
            stack_order.append(order_count)
            # bp_index += 1
            recursive_seq(hash, stack_component, stack_order, target, const_new_order_count, bp_index + 1, cnt)
    else:
        if bp_index < 1:
            return 0
        else:
            recursive_seq(hash, stack_component, stack_order, stack_component[bp_index-1], stack_order[bp_index-1] +1, bp_index - 1, cnt)
        
    return size(stack_component)



# sc = []; so = []
# tmp_hash = loadtxt('tmp.hash')
# print 0
# size = recursive_seq(hash=tmp_hash, stack_component=sc, stack_order=so, index=0)
