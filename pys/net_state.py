#!/bin/python
# this script aims to draw rasterogram of neural network
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
parser.add_argument('dir', type = str, help = 'directory of source data and output data')
args = parser.parse_args()

# Run simulations

#subprocess.call(['cp', 'doc/config_net.ini', args.dir])
#p = subprocess.call(['./bin/net.out', 'doc/config_net.ini', args.dir])
p = 0
if p == 0:
    typepath = args.dir + '/neuron_type.csv'
    types = np.genfromtxt(typepath, delimiter = ',', dtype = 'b')
    types = types[:-1]

    # import config file

    config = cp.ConfigParser()
    config.read(args.dir + '/config_net.ini')
    neuron_num = int(config.get('Network Parameters', 'NeuronNumber'))
    tmax = float(config.get('Time', 'MaximumTime'))

    # create canvas

    fig = plt.figure(figsize = (14,4), dpi = 60)

    # config the plotting range (unit millisecond)
    dt = 5
    t_start = tmax - 1500 
    t_end = tmax 
    t = np.arange(t_start, t_end, dt)
    counts_exc = np.zeros(len(t))
    counts_inh = np.zeros(len(t))

    # draw raster plot

    ax1 = plt.subplot2grid((2,6), (0,0), colspan = 4, rowspan = 1)
    f = open(args.dir + '/raster.csv')
    counter = 1
    mrate = np.zeros(neuron_num)
    isi_e = np.array([])
    isi_i = np.array([])
    for line in f:
        spike_str = line.replace('\n', '').split(',')
        spike_str.remove('')
        spikes = [float(i) for i in spike_str if float(i) < tmax]
        mrate[counter-1] = len(spikes)
        counter_list = np.ones(len(spikes)) * counter
        if types[counter - 1]:
            ax1.scatter(spikes, counter_list, s = 1, c = 'b')
            isi_e = np.append(isi_e, np.diff(spikes))
        else:
            ax1.scatter(spikes, counter_list, s = 1, c = 'r')
            isi_i = np.append(isi_i, np.diff(spikes))
        counter += 1
    f.close();
    ax1.set_xlim(t_start, t_end)
    ax1.set_ylim(0, counter + 1)
    ax1.set_ylabel('Indices')
    ax1.set_xlabel('Time (ms)')
    ax1.grid(linestyle='--')

    # plot the histogram of mean firing rate

    ax2 = plt.subplot2grid((2,6), (1,0), colspan = 2, rowspan = 1)
    ax2.hist(mrate[(types==1).nonzero()]/(mrate.sum()/neuron_num), 50, label = 'EXC Neurons')
    ax2.hist(mrate[(types==0).nonzero()]/(mrate.sum()/neuron_num), 50, label = 'INH Neurons')
    #ax2.hist(mrate/(mrate.sum()/neuron_num), 50)
    ax2.set_xlabel('Rate/mean rate')
    ax2.set_ylabel('Density')
    ax2.grid(linestyle='--')
    ax2.legend()
 
    # draw mean firing rate of network

    ax3 = plt.subplot2grid((2,6), (1,2), colspan = 2, rowspan = 1)
    ax3.hist(isi_e, 100, label = 'EXC Neurons')
    ax3.hist(isi_i, 100, label = 'INH Neurons')
    ax3.set_xlabel('ISI (ms)')
    ax3.set_ylabel('Density')
    ax3.legend()
    ax3.grid(linestyle='--')

    # calculate the mutual information between neurons and LFP
    X,Y = np.meshgrid(range(neuron_num + 1),range(neuron_num + 1))
    mat = np.genfromtxt(args.dir + '/mat.csv', delimiter = ',', dtype = int)
    mat = mat[:,:-1]

    ax4 = plt.subplot2grid((2,6), (0,4), colspan = 2, rowspan = 2)
    cax4 = ax4.pcolormesh(X, Y, mat, cmap = 'Greys')
    ax4.set_xlabel('Post-neurons')
    ax4.set_ylabel('Pre-neurons')
    ax4.set_title('Connectivity Mat')
    #ax4.grid(linestyle = '--')
    fig.colorbar(cax4, ax = ax4)

    # tight up layout and save the figure

    plt.tight_layout()
    plt.savefig(args.dir + '/net_state.png')
    plt.close()
