#!/usr/bin/python
"""
This script aims to calculate the autocovariance of long continues data with built-in fft functions;
    auto-covariance = np.fft.fft( np.fft.ifft(data) * np.conj( np.fft.ifft(data) ) )
"""
import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import struct

f = open('./data/tmp/singleI.bin', 'rb')
sshape = f.read(2*8)
shape = struct.unpack('Q'*2, sshape)
sdata = f.read(int(shape[1])*8)
f.close()
data = np.array(struct.unpack('d'*int(shape[1]), sdata))

# cov = np.correlate(data,data,'same')
fdata = np.fft.ifft(data)
fdata = fdata*np.conj(fdata)
cov = np.fft.fft(fdata)
# print cov
pltcov = abs(cov[0:8000]) - data.mean()**2
pltcov = pltcov/pltcov.max()
plt.plot(np.arange(len(pltcov))*0.03125, pltcov)
plt.xlabel('Time(ms)')
plt.ylabel('Normalized Autocovariance')
plt.grid()
plt.show()
