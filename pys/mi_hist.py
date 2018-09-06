# -*- coding: utf-8 -*-
#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import sys
import argparse

parser = argparse.ArgumentParser(description = "Plot the histogram of maximum mutual information based on pairwised analysis.")
parser.add_argument('filename', metavar = 'target_data', type = str, nargs = 1, help = 'path of mi data')
parser.add_argument('mat', metavar = 'target_con_mat', default = 'mat.csv', type = str, nargs = 1, help = 'path of connectivity matrix')
parser.add_argument('-s', '--show', action = 'store_true', help = 'show the plotted figure')
args = parser.parse_args()
arr = np.genfromtxt(args.filename[0], delimiter = ',')
mat = np.genfromtxt(args.mat[0], delimiter = ',')
mat = np.delete(mat, -1, 1)
data_len = 1e8
# examine the size of data matrix
(nrow, ncol) = arr.shape

one = np.array([])
negone = np.array([])
zero = np.array([])
for i in range(nrow):
    for j in range(ncol):
        if (i != j):
            if mat[i][j] == 1:
                one = np.append(one, arr[i][j])
            else:
                if mat[j][i] == 1:
                    negone = np.append(negone, arr[i][j])
                else:
                    zero = np.append(zero, arr[i][j])

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
bin_num = 60
one = np.array([np.log10(ele*1.0) for ele in one if ele != 0])
negone = np.array([np.log10(ele*1.0) for ele in negone if ele != 0])
zero = np.array([np.log10(ele*1.0) for ele in zero if ele != 0])
one_counts, one_edge = np.histogram(one, int(bin_num*0.3))
negone_counts, negone_edge = np.histogram(negone, int(bin_num*0.4))
zero_counts, zero_edge = np.histogram(zero, int(bin_num*0.8))
one_counts = one_counts * 1.0 / (nrow*(ncol-1)) / (one_edge[1] - one_edge[0])
negone_counts = negone_counts * 1.0 / (nrow*(ncol-1)) / (negone_edge[1] - negone_edge[0])
zero_counts = zero_counts * 1.0 / (nrow*(ncol-1)) / (zero_edge[1] - zero_edge[0])

#one_edge = 10**one_edge
#zero_edge = 10**zero_edge
#plt.semilogx(one_edge, one_counts, label = 'connected')
#plt.semilogx(zero_edge, zero_counts, label = 'disconneted')
#plt.bar(one_edge[:-1], one_counts, width = one_edge[1] - one_edge[0], label = 'connected')
#plt.bar(negone_edge[:-1], negone_counts, width = negone_edge[1] - negone_edge[0], label = 'inversely connected')
#plt.bar(zero_edge[:-1], zero_counts, width = zero_edge[1] - zero_edge[0], label = 'disconnected')
one_edge = (one_edge[:-1] + one_edge[1:])/2
negone_edge = (negone_edge[:-1] + negone_edge[1:])/2
zero_edge = (zero_edge[:-1] + zero_edge[1:])/2
plt.figure(figsize = (14,10), dpi = 80)
plt.plot(one_edge, one_counts, '.-', label = 'connected')
plt.plot(negone_edge, negone_counts, '.-', label = 'inversely connected')
plt.plot(zero_edge, zero_counts, '.-', label = 'disconnected')
plt.axvline(-np.log10(data_len*1.0), color='cyan')

# plot settings:
xmajorLocator   = MultipleLocator(1) 
xmajorFormatter = FormatStrFormatter('$10^{%d}$') # set formats of x axis
xminorLocator   = MultipleLocator(0.2) 

ymajorLocator   = MultipleLocator(1e-2) 
#ymajorFormatter = FormatStrFormatter('%1.1f') # set formats of y axis
yminorLocator   = MultipleLocator(0.2e-2) 

plt.xlabel('Mutual Information', fontsize = 20)
plt.ylabel('Probability Density', fontsize = 20)
ax = plt.gca()
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)

ax.yaxis.get_major_formatter().set_powerlimits((0,1))
#ax.yaxis.set_major_locator(ymajorLocator)
#ax.yaxis.set_minor_locator(yminorLocator)

plt.legend(fontsize = 20)
plt.grid()
if args.show:
    plt.show()
else:
    plt.savefig('mi_hist.png')
    plt.close()
