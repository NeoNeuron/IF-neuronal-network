#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

data = np.loadtxt('./mi_max.csv', delimiter = ',')
plt.figure(figsize = (12,18), dpi = 80, facecolor = None)
plt.subplot(3,1,1)
plt.imshow(data, aspect = 'auto', extent = [10000, 590000, 100, 5], cmap = cm.jet)
plt.xlabel('Max time(ms)', fontsize = 18)
plt.ylabel('No.of bins in histograms', fontsize = 18)
plt.title('Maximum mutual information')
plt.colorbar()
data = np.loadtxt('./snr.csv', delimiter = ',')
plt.subplot(3,1,2)
plt.imshow(data, aspect = 'auto', extent = [10000, 590000, 100, 5], cmap = cm.jet)
plt.xlabel('Max time(ms)', fontsize = 18)
plt.ylabel('No.of bins in histograms', fontsize = 18)
plt.title('Signal-to-Noise Ratio')
plt.colorbar()
data = np.loadtxt('./peak_dev.csv', delimiter = ',')
plt.subplot(3,1,3)
plt.imshow(data, aspect = 'auto', extent = [10000, 590000, 100, 5], cmap = cm.jet)
plt.xlabel('Max time(ms)', fontsize = 18)
plt.ylabel('No.of bins in histograms', fontsize = 18)
plt.title('Deviation of maximum mutual information')
plt.colorbar()
plt.savefig('./fig.png')
plt.close()
