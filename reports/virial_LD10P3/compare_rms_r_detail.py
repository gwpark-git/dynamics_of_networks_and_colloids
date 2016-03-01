
NP_arr = asarray([400, 500, 600, 700])

plt.clf()
plt.ion()
plt.grid()
plt.plot(NP_arr, mean_Rrms_arr, 'bo-', label = 'average RMS for distance of bridges')
plt.xticks(NP_arr)
plt.yticks(mean_Rrms_arr, ['%.3f'%val for val in mean_Rrms_arr])
plt.xlabel('number of particls')
plt.ylabel('average RMS for distance of bridges')
plt.legend(loc = 'lower right', numpoints=1)
plt.axis([350, 750, 0.973, 0.985])
plt.show()

