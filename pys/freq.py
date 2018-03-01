#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct

f = open('./data/tmp/singleI.bin', 'rb')
sshape = f.read(2*8)
shape = struct.unpack('Q'*2, sshape)
sdata = f.read(shape[1]*8)
f.close()
data = np.array(struct.unpack('d'*shape[1], sdata))
# calculate spectrum;
fdata = np.fft.fftshift(np.fft.fft(data))
fdata = abs(fdata)

plt.loglog(fdata[int(len(fdata)/2):-1])
plt.show()
