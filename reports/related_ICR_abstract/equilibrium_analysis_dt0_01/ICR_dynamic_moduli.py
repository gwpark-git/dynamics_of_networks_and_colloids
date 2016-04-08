
# def get_dynamic_moduli(t, Gt):
#     w = sort(1./t)
#     Gp = zeros(size(w))
#     Gpp = zeros(size(w))
#     for i in range(size(w)):
#         if i%1000==0:
#             print i, w[i]
#         for j in range(1,size(t)):
#             Gp[i] += 0.5*(t[j] - t[j-1])*(Gt[j] + Gt[j-1])*sin(w[i]*t[j])
#             Gpp[i] += 0.5*(t[j] - t[j-1])*(Gt[j] + Gt[j-1])*cos(w[i]*t[j])
#         Gp[i] *= w[i]
#         Gpp[i] *= w[i]
#     return w, Gp, Gpp

w_NP0400, Gp_NP0400, Gpp_NP0400 = get_dynamic_moduli(t_NP0400[:100], acf_NP0400[:100])
w_NP1000, Gp_NP1000, Gpp_NP1000 = get_dynamic_moduli(t_NP1000[:300], acf_NP1000[:300])




plt.close()
plt.ion()
plt.loglog(w_NP0400, Gp_NP0400, 'b-')
plt.loglog(w_NP0400, Gpp_NP0400, 'r-')

plt.loglog(w_NP1000, Gp_NP1000, 'b--')
plt.loglog(w_NP1000, Gpp_NP1000, 'r--')

plt.grid()
plt.show()
