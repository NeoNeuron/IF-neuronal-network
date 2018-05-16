#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import struct as st
import sys

path = sys.argv[1]
V_ind = int(sys.argv[2])
I_ind = int(sys.argv[3])
raster_ind = int(sys.argv[4])
dt = 0.5
dat_len = 4000
# load Voltage and Current
f = open(path + 'V.bin', 'rb')
V_shape = st.unpack('Q'*2, f.read(8*2))
V = np.array(st.unpack('d'*dat_len*V_shape[1], f.read(8*dat_len*V_shape[1])))
f.close()
V = V.reshape((dat_len,V_shape[1]))
V_single = V[:,V_ind]

f = open(path + 'I.bin', 'rb')
I_shape = st.unpack('Q'*2, f.read(8*2))
I = np.array(st.unpack('d'*dat_len*I_shape[1], f.read(8*dat_len*I_shape[1])))
f.close()
I = I.reshape((dat_len,I_shape[1]))
I_single = I[:,I_ind]

# load raster data:
raster = np.genfromtxt(path + 'raster.csv', delimiter = ',', skip_header = raster_ind, skip_footer = V_shape[1] - raster_ind - 1)

# plot v and I
plt.figure(figsize=(20,10), dpi = 80)
plt.subplot(2,1,1)
plt.plot((np.arange(dat_len)+1)*dt, V_single, label = 'V')
plt.xlabel('time(ms)')
plt.ylabel('Potential')
plt.legend()
for i in raster:
    if i <= dat_len*dt:
        plt.axvline(i, color='red')
    else:
        break
plt.subplot(2,1,2)
plt.plot((np.arange(dat_len)+1)*dt, I_single, label = 'I')
plt.xlabel('time(ms)')
plt.ylabel('Current')
plt.legend()
for i in raster:
    if i <= dat_len*dt:
        plt.axvline(i, color='red')
    else:
        break
# savefig
plt.savefig('V_I.png')
#plt.show()
plt.close()
