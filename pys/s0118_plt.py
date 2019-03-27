# this script aims to draw the compact figure based on s0118.py with existing data
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import argparse
import configparser as cp
import struct as st
import subprocess

# config input parameter:

parser = argparse.ArgumentParser(description = "Integrated script of to anaylsis spike-LFP mutual information in 2-D neuronal net.")
parser.add_argument('dir', type = str, help = 'directory of source data')
args = parser.parse_args()

typepath = args.dir + '/neuron_type.csv'
types = np.genfromtxt(typepath, delimiter = ',', dtype = 'b')
types = types[:-1]

# import config file

config = cp.ConfigParser()
config.read('doc/config_net.ini')
neuron_num = int(config.get('Network Parameters', 'NeuronNumber'))
tmax = float(config.get('Time', 'MaximumTime'))

# create canvas

fig = plt.figure(figsize = (12,10), dpi = 60)

# config the plotting range (unit millisecond)
dt = 5
t_start = 2500
t_end = 3000
t = np.arange(t_start, t_end, dt)
counts_exc = np.zeros(len(t))
counts_inh = np.zeros(len(t))

# draw raster plot

ax1 = plt.subplot2grid((12,5), (0,0), colspan = 3, rowspan = 3)
f = open(args.dir + '/raster.csv')
counter = 1
for line in f:
    spike_str = line.replace('\n', '').split(',')
    spike_str.remove('')
    spikes = [float(i) for i in spike_str if float(i) < tmax]
    counter_list = np.ones(len(spikes)) * counter
    if types[counter - 1]:
        ax1.scatter(spikes, counter_list, s = 1, c = 'b')
        for ele in spikes:
            if ele >= t_start and ele < t_end:
                counts_exc[int((ele - t_start)/dt)] += 1
            elif ele >= t_end:
                break
    else:
        ax1.scatter(spikes, counter_list, s = 1, c = 'r')
        for ele in spikes:
            if ele >= t_start and ele < t_end:
                counts_inh[int((ele - t_start)/dt)] += 1
            elif ele >= t_end:
                break
    counter += 1
f.close();
ax1.set_xlim(t_start, t_end)
ax1.set_ylim(0, counter + 1)
ax1.set_ylabel('Indices')
ax1.grid(linestyle='--')

# draw mean firing rate of network

ax2 = plt.subplot2grid((12,5), (3,0), colspan = 3, rowspan = 3)
ax2.plot(t + dt/2, counts_exc, linestyle = '--', linewidth = 1, label='exc neuron')
ax2.plot(t + dt/2, counts_inh, linestyle = '--', linewidth = 1, label='inh neuron')
ax2.plot(t + dt/2, counts_exc + counts_inh, linewidth = 1, label='tot neuron')
ax2.set_ylabel('Mean Firing Counts')
ax2.legend()
ax2.grid(linestyle='--')

# draw LFP trajectory

config_lfp = cp.ConfigParser()
config_lfp.read('doc/config_lfp.ini')
lfp_tmin = float(config_lfp.get('Time', 'TimeRangeMin'))

f = open(args.dir + '/lfp.bin', 'rb')
shape = st.unpack('QQ', f.read(16))
lfp = np.array(st.unpack('d'*shape[1], f.read(8*shape[1])))
f.close()
# slicing the original lfp series;
dt_sampling = 0.5
ax3 = plt.subplot2grid((12,5), (6,0), colspan = 3, rowspan = 3)
ax3.plot(np.arange(t_start, t_end, dt_sampling), lfp[int((t_start-lfp_tmin)/dt_sampling):int((t_end-lfp_tmin)/dt_sampling)], label = "LFP")
ax3.set_xlabel('Time(ms)')
ax3.set_ylabel('LFP')
ax3.legend()
ax3.grid(linestyle = '--')

# draw power spectrum of LFP

lfp_fft = abs(np.fft.fft(lfp))

# config frequency range (unit Hz)

f_start = 0
f_end = 100
df = 1e3/tmax
freq = np.arange(f_start, f_end, df)

ax4 = plt.subplot2grid((12,5), (9,0), colspan = 3, rowspan = 3)
ax4.semilogy(freq, lfp_fft[0 : 0 + len(freq)], linewidth = 0.1)
ax4.set_xlabel('Freq(Hz)')
ax4.set_ylabel('LFP Power Spectra')
ax4.grid(linestyle='--')

# import mutual information data

mi_dat = np.genfromtxt(args.dir + '/mi_lfp_scan.csv', delimiter = ',')

# draw grid scatter plot of maximum mutual information

xn_gd = 10
yn_gd = 10
X,Y = np.meshgrid(range(xn_gd + 1),range(yn_gd + 1))
X = X / xn_gd
Y = Y / yn_gd
ax5 = plt.subplot2grid((12,5), (0,3), colspan = 2, rowspan = 4)
cax5 = ax5.pcolormesh(X, Y, np.reshape(mi_dat.max(1), (xn_gd,yn_gd)), cmap = 'jet', norm = colors.LogNorm(vmin = mi_dat.max(1).min(), vmax = mi_dat.max(1).max()))
ax5.set_xlabel('X')
ax5.set_ylabel('Y')
ax5.set_title('Max MI')
ax5.grid(linestyle = '--')
fig.colorbar(cax5, ax = ax5)

# draw the position of electrode
xe = float(config_lfp.get('electrode', 'PosX'))
ye = float(config_lfp.get('electrode', 'PosY'))
ax5.scatter(xe, ye, s = 60, marker = "+", c = 'y')

# draw grid scatter plot of maximum time delay

ax6 = plt.subplot2grid((12,5), (4,3), colspan = 2, rowspan = 4)
mi_delay = (mi_dat.argmax(1)-15)*0.5
if abs(mi_delay.max()) > abs(mi_delay.min()):
    cmax = abs(mi_delay.max())
    cmin = -cmax
else:
    cmax = abs(mi_delay.min())
    cmin = -cmax
cax6 = ax6.pcolormesh(X, Y, np.reshape(mi_delay, (xn_gd,yn_gd)), cmap = 'seismic', norm = colors.Normalize(vmin = cmin, vmax = cmax))
ax6.set_xlabel('X')
ax6.set_ylabel('Y')
ax6.set_title('Max Time Delay')
ax6.grid(linestyle = '--')
cb6 = fig.colorbar(cax6, ax = ax6)
cb6.set_label('delay time (ms)')
ax6.scatter(xe, ye, marker = "+", c = 'y')

# draw grid scatter plot of neuronal type

ax7 = plt.subplot2grid((12,5), (8,3), colspan = 2, rowspan = 4)

# create colorbar:

cm7 = colors.ListedColormap(['red', 'blue'])
norm7 = colors.BoundaryNorm([0,0.5,1], cm7.N)
cax7 = ax7.pcolormesh(X, Y, np.reshape(types, (xn_gd,yn_gd)), cmap = cm7, norm = norm7)
ax7.set_xlabel('X')
ax7.set_ylabel('Y')
ax7.set_title('Neuronal Types')
ax7.grid(linestyle = '--')
cb7 = fig.colorbar(cax7, ax = ax7)
cb7.set_label('neuronal type')

# tight up layout and save the figure

plt.tight_layout()
plt.savefig(args.dir + '/raster.eps')
plt.close()
