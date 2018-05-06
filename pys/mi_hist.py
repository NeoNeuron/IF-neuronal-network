#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import sys
arr_name = sys.argv[1]
arr = np.genfromtxt(arr_name, delimiter = ',')
mat = np.genfromtxt('mat.csv', delimiter = ',')
mat = np.delete(mat, -1, 1)

one = np.zeros(10000)
zero = np.zeros(10000)
one_counter = 0
zero_counter = 0
for i in range(100):
    for j in range(100):
        if (i != j):
            if mat[i][j] == 1:
                one[one_counter] = arr[i][j]
                one_counter += 1
            else:
                zero[zero_counter] = arr[i][j]
                zero_counter += 1
one = one[:one_counter-1]
zero = zero[:zero_counter-1]
# linearly binned:
#bin_num = 100
##binsize = 1e-7
##one_bin = int(np.ceil((one.max() - one.min()) / binsize))
##zero_bin = int(np.ceil((zero.max() - zero.min()) / binsize))
#one_counts, one_edge = np.histogram(one, bin_num)
#zero_counts, zero_edge = np.histogram(zero, bin_num)
#one_counts = one_counts * 1.0 / sum(one_counts)
#zero_counts = zero_counts * 1.0 / sum(zero_counts)
#one_edge = (one_edge[:-1] + one_edge[1:])/2
#zero_edge = (zero_edge[:-1] + zero_edge[1:])/2
# loganormially binned:
bin_num = 100
one = np.array([np.log10(ele*1.0) for ele in one if ele != 0])
zero = np.array([np.log10(ele*1.0) for ele in zero if ele != 0])
one_counts, one_edge = np.histogram(one, bin_num)
zero_counts, zero_edge = np.histogram(zero, bin_num)
one_counts = one_counts * 1.0 / sum(one_counts)
zero_counts = zero_counts * 1.0 / sum(zero_counts)
#one_edge = (one_edge[:-1] + one_edge[1:])/2
#zero_edge = (zero_edge[:-1] + zero_edge[1:])/2
#one_edge = 10**one_edge
#zero_edge = 10**zero_edge
#plt.semilogx(one_edge, one_counts, label = 'connected')
#plt.semilogx(zero_edge, zero_counts, label = 'disconneted')
plt.bar(one_edge[:-1], one_counts, width = one_edge[1] - one_edge[0], label = 'connected')
plt.bar(zero_edge[:-1], zero_counts, width = zero_edge[1] - zero_edge[0], label = 'disconnected')
# plot settings:
xmajorLocator   = MultipleLocator(1) 
xmajorFormatter = FormatStrFormatter('$10^{%1.1f}$') # set formats of x axis
xminorLocator   = MultipleLocator(0.2) 

ymajorLocator   = MultipleLocator(1e-2) 
#ymajorFormatter = FormatStrFormatter('%1.1f') # set formats of y axis
yminorLocator   = MultipleLocator(0.2e-2) 

plt.xlabel('Mutual information')
plt.ylabel('Probability')
ax = plt.gca()
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)

ax.yaxis.get_major_formatter().set_powerlimits((0,1))
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_minor_locator(yminorLocator)

plt.legend()
plt.grid()
if len(sys.argv) == 2:
    plt.show()
elif len(sys.argv) == 3:
    plt.savefig(sys.argv[2])
    plt.close()
