#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

data = np.loadtxt('./mi_max.csv', delimiter = ',')
plt.figure(figsize = (12,8), dpi = 80, facecolor = None)
# plt.subplot(3,1,1)
im = plt.imshow(data, aspect = 'auto', extent = [30000, 300001, 20, 3], cmap = cm.jet)
plt.xlabel('Max time(ms)', fontsize = 18)
plt.ylabel('No.of bins in histograms', fontsize = 18)
bfmt = mpl.ticker.ScalarFormatter()
bfmt.set_powerlimits((0,1))
bar = plt.colorbar(im, format = bfmt)
bar.set_label('Maximum Mutual Information')
# data = np.loadtxt('./snr.csv', delimiter = ',')
# plt.subplot(3,1,2)
# plt.imshow(data, aspect = 'auto', extent = [30000, 300001, 20, 3], cmap = cm.jet)
# plt.xlabel('Max time(ms)', fontsize = 18)
# plt.ylabel('No.of bins in histograms', fontsize = 18)
# plt.title('Signal-to-Noise Ratio')
# plt.colorbar()
# data = np.loadtxt('./peak_dev.csv', delimiter = ',')
# plt.subplot(3,1,3)
# plt.imshow(data, aspect = 'auto', extent = [30000, 300001, 20, 3], cmap = cm.jet)
# plt.xlabel('Max time(ms)', fontsize = 18)
# plt.ylabel('No.of bins in histograms', fontsize = 18)
# plt.title('Deviation of maximum mutual information')
# plt.colorbar()
plt.savefig('fig.png')
plt.close()
