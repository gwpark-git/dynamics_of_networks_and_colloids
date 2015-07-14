# this package is curve-fitting using Chebyshev polynomial 1st kind and then generating polynomial coefficient
# The methodology is developed by Prof. Kwang Soo Cho
# The source code is developed by Gun Woo Park


from numpy import *
from scipy import linalg



def factorial(N):
    re = 1
    for i in range(1,N+1):
        re *= i
    return re

def Chebyshev_1st(n, x):
    if n==0:
        return 1.0
    elif n==1:
        return x
    return 2.0*x*Chebyshev_1st(n-1, x) - Chebyshev_1st(n-2, x)

def coe_Chebyshev_1st(k, n):
    if k<0 or k > n or (k==0 and n==1):
        return 0
    elif (k==1 and n==1) or (k==0 and n==0):
        return 1
    else:
        return 2*coe_Chebyshev_1st(k-1, n-1) - coe_Chebyshev_1st(k, n-2)

def cal_Cheby_1st(cn, x):
    N_n = size(cn)
    N_x = size(x)
    y = zeros(N_x)
    for i in range(N_x):
        for k in range(N_k):
            y[i] += Chebyshev_1st(k, x)*x
    return 0

# def L_nk(N, t):
#     L = zeros([N, N])
#     for i in range(N):
#         for j in range(N):
#             for k in range(size(t)):
#                 L[i,j] += Chebyshev_1st(i, t[k])*Chebyshev_1st(j, t[k])
#     return L_nk


def get_Cheby_1st_coe(t, y, N):
    a = zeros(N) # coefficient
    L = zeros([N, N]) #La = yT
    for i in range(N):
        for j in range(N):
            for k in range(size(t)):
                L[i,j] += Chebyshev_1st(i, t[k])*Chebyshev_1st(j, t[k])

    yT = zeros(N)
    for i in range(N):
        for j in range(size(y)):
            yT[i] += y[j]*Chebyshev_1st(i, t[j])
    a = linalg.inv(L).dot(yT)
    return a

def get_polynomial_from_Cheby(x, y, N):
    # note that dx must be constant for this method
    x_min = min(x)
    x_max = max(x)
    x_c = (max(x) - min(x))/2.
    dx = x[1] - x[0]
    t = 2.*(x - x_c)/dx
    a = get_Cheby_1st_coe(t, y, N)
    
    b = zeros(N)

    for n in range(N):
        for k in range(n, N):

            b[n] += coe_Chebyshev_1st(k, n)*a[k]
    c = zeros(N)
    for n in range(N):
        for k in range(n, N):
            c[n] += ((2./dx)**k)*(factorial(k)/(factorial(n)*factorial(k-n)))*((-x_c)**(k-n))*b[k]
    return c

    
def gen_from_polynomial(x, c):
    N_x = size(x)
    y = zeros(N_x)
    for i in range(N_x):
        for k in range(size(c)):
            y[i] += c[k]*x[i]**k
    return y

def diff_gen_from_polynomial(x, c):
    N_x = size(x)
    y = zeros(N_x)
    for i in range(N_x):
        for k in range(size(c)-1):
            y[i] += c[k+1]*float(k+1)*x[i]**(k)
    return y
