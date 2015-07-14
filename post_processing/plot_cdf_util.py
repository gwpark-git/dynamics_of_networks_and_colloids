
# c = get_polynomial_from_Cheby(tnew, ynew, 8)
# t_cheby = linspace(min(tnew), max(tnew), 100)
# y_cheby = gen_from_polynomial(t_cheby, c)
# diff_cheby = diff_gen_from_polynomial(t_cheby, c)
# def U_ij(r):
#     return (1./3.)*(1. - r)**2.*(2. + r)

ref_line = asarray([[1., -0.1],
                    [1., 5]])


plt.clf()
fig, ax1 = plt.subplots()
leg = []
ax2 = ax1.twinx()
leg.append(ax1.plot(t, y, 'bo', markerfacecolor='white', markeredgecolor='blue', label='CDF(r), data')[0])
leg.append(ax1.plot(tnew, ynew, 'r-', label='CDF(r), interpolated')[0])
# leg.append(ax1.plot(t_cheby, y_cheby, 'g.', label='CDF(r), Chebyshev')[0])
leg.append(ax2.plot(tnew, diff, 'k-', label = 'PDF(r)')[0])
# leg.append(ax2.plot(tnew, 4.*exp(-U_ij(tnew)), 'g-', label = 'U(r)')[0])

# leg.append(ax2.plot(t_cheby, diff_cheby, 'g-', label='PDF(r), Chebyshev')[0])
leg.append(ax2.plot(ref_line[:,0], ref_line[:,1], 'r--', linewidth=3, label='Range for Repulsion')[0])

ax1.grid('on')
ax1.set_xlabel('distance, r')
ax1.set_ylabel('cumulative distribution function up to r=1.5, CDF(r)')
ax2.set_ylabel('probability distribution function from cdf, PDF(r)')
# legend_line = ax.legend(handles=leg, loc = 'upper left', numpoints=1)
plt.legend(handles=leg, loc = 'upper left', numpoints=1, ncol=2)
ax1.axis([t_min, t_max, -0.1, 1])
ax2.axis([t_min, t_max, -0.1, 5])
plt.savefig('Np100_repulsion.pdf')
plt.show()
