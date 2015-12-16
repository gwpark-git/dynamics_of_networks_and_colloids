
from itertools import combinations
from itertools import combinations_with_replacement

from numpy import *

def get_combinations_with_replacement(Np, Nc):
    return asarray(tuple(combinations_with_replacement(arange(Np), Nc)))

def get_weight_combinations(Np, Nc):
    tmp_w = get_combinations_with_replacement(Np, Nc)
    Ns = shape(tmp_w)[0]
    w = zeros([Ns, Np])
    for s in range(Ns):
        for i in range(Nc):
            w[s, tmp_w[s, i]] += 1
    return w
