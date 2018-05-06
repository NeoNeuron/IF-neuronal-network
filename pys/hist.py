#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct as st
import sys
#import data:
f = open(sys.argv[1], 'rb')
size = [0,0]
size[0] = st.unpack('Q', f.read(8))[0]
size[1] = st.unpack('Q', f.read(8))[0]
#prepare data container
dat = np.empty(size[1])
for i in range(size[1]):
    dat[i] = st.unpack('d', f.read(8))[0]
f.close()
print dat.min()
#calculate histogram
bins = int((dat.max() - dat.min())/float(sys.argv[2]))
counts, edges = np.histogram(dat, bins = bins)
counts = counts * 1.0 / sum(counts)
width = edges[1] - edges[0]
# logscale
if sys.argv[3] == 'log':
    for i in range(len(counts)):
        if counts[i] != 0:
            counts[i] = np.log10(counts[i]*1.0)
    plt.bar(edges[:-1], counts, width = width)
    plt.xlabel('Value of current')
    plt.ylabel('Counts (log10-scale)')
    plt.grid()
elif sys.argv[3] == 'lin':
    plt.bar(edges[:-1], counts, width = width)
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
    plt.xlabel('Value of current')
    plt.ylabel('Counts (linear-scale)')
    plt.grid()
# adjust fig ouput:
if len(sys.argv) == 4:
    plt.show()
elif len(sys.argv) == 5:
    plt.savefig(sys.argv[4])
    plt.close()
