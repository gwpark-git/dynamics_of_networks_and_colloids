
NP_arr = asarray([400, 500, 600, 700])

plt.clf()
plt.ion()
plt.grid()
plt.plot(NP_arr, mean_beta_arr, 'bo-', label = 'average detachment frequency in equilibrium')
plt.xticks(NP_arr)
plt.yticks(mean_beta_arr, ['%.3f'%val for val in mean_beta_arr])
plt.xlabel('number of particls')
plt.ylabel('average detachment frequency')
plt.legend(loc = 'lower right', numpoints=1)
plt.axis([350, 750, 1.02, 1.06])
plt.show()

