import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
dt = 0.5
# Calculating Mutual information:
ntd = 15
ptd = 15
drange = str(ntd) + ',' + str(ptd)
path = sys.argv[1]
spike_ind = sys.argv[2]
lfp_ind = sys.argv[3]
binsize = sys.argv[4]
trange = '5e2,1e6'
program = './bin/mi_bd_unity.out'
subprocess.call([program, path, spike_ind, lfp_ind, trange, str(dt), drange, binsize])
data = np.genfromtxt(path + 'mi_bd.csv', delimiter = ',')
# Plotting:
plotway = 'linear'
if plotway == 'linear':
    plt.plot(data[:,0]*dt, data[:,1], '-*', label = 'mi')
    plt.plot(data[:,0]*dt, data[:,2], '-*', label = 'mi_shuffle')
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
    plt.semilogy(data[:,0]*dt, data[:,1], '-*', label = 'mi')
    plt.semilogy(data[:,0]*dt, data[:,2], '-*', label = 'mi_shuffle')

plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.legend()
plt.grid()
plt.savefig('mi')
plt.close()
