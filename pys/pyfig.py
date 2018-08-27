#!/usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

arr = np.genfromtxt('data/dataRepo/1.3kHz/raster.csv', delimiter = ',', skip_footer = 1)
arr = np.delete(arr, -1, 0)
# Set the time step and time range;
dt = 0.5
tmax = 900000000
spike_cum = np.empty(int(tmax/dt) + 1)
spike_cum[0] = 0
counter = 0
print arr[-1]
print arr.max()
print len(np.arange(0, tmax + dt, dt))
print len(spike_cum)
for i in range(len(spike_cum) - 1):
  if counter < len(arr):
    if dt*(i + 1) >= arr[counter]:
      spike_cum[i + 1] = spike_cum[i] + 1
      counter += 1
    else:
      spike_cum[i + 1] = spike_cum[i]
  else:
    spike_cum[i + 1] = spike_cum[i]

print counter
print len(arr)
plt.plot(np.arange(0, tmax + dt, dt), spike_cum)
plt.xlabel('Time')
plt.ylabel('Spike number')
plt.grid()
plt.savefig('fig.png')
plt.close()
