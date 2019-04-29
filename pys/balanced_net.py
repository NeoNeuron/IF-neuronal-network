#!/bin/python
# draw mean firing rate as a funtion of network inputs
import numpy as np
import matplotlib.pyplot as plt
import argparse
import configparser as cp
import struct as st
import subprocess

# rewrite the configparser.ConfigParser
class MyConfigParser(cp.ConfigParser):
    def __init__(self,defaults=None):
        cp.ConfigParser.__init__(self,defaults=None)
    def optionxform(self, optionstr):
        return optionstr
#---------------------
# config input parameter:

parser = argparse.ArgumentParser(description = "Integrated script for study the input-output relation of balanced network.")
parser.add_argument('dir', type = str, help = 'directory of source data and output data')

args = parser.parse_args()

config = MyConfigParser()
config.read('doc/config_net.ini')
config['Output']['SaveV'] = 'false'
config['Output']['SaveI'] = 'false'
config['Output']['SaveGE'] = 'false'
config['Output']['SaveGI'] = 'false'
config['Time']['MaximumTime'] = '1e3'

# iteration range:

pr = float(config['Driving Settings']['pr'])
ps = np.arange(1e-4, 5.1e-3, 4e-4)
frate = np.zeros(len(ps))
for i in range(len(ps)):
    config['Driving Settings']['ps'] = str(ps[i])
    # save config file
    with open(args.dir + '/config_net.ini', 'w') as configfile:
        config.write(configfile)
    p = subprocess.call(['./bin/net.out', args.dir+ '/config_net.ini', args.dir])
    if p == 0:
        f = open(args.dir + '/raster.csv')
        for line in f:
            frate[i] += len(line.split(',')) - 1

frate /= 1e-3 * float(config['Time']['MaximumTime']) * int(config['Network Parameters']['NeuronNumber'])

# plot figure
fig = plt.figure(figsize=(8,6), dpi =96)
ax = fig.subplots(1,1)
ax.plot(pr*ps, frate, '-o')
ax.set_xlabel(r'$f\nu(ms^{-2})$', fontsize = 14)
ax.set_ylabel('Mean firing rate(Hz)', fontsize = 14)
ax.grid(linestyle = '--')
ax.xaxis.get_major_formatter().set_powerlimits((0,1))
plt.tight_layout()
plt.savefig('bnet.png')
plt.show()
