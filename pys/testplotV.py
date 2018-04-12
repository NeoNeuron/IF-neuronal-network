#!/usr/bin/python
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import struct

V = np.empty(5000)
I = np.empty(5000)
f = open('data/dataRepo/net1.3kHz/V.bin', 'rb')
f.seek(16 + 2000000000*8)
for i in range(len(V)):
    mid = f.read(8)
    V[i] = struct.unpack('d', mid)[0]
    f.seek(8,1)
f.close()

f = open('data/dataRepo/net1.3kHz/I.bin', 'rb')
f.seek(16 + 2000000000*8+8)
for i in range(len(I)):
    mid = f.read(8)
    I[i] = struct.unpack('d', mid)[0]
    f.seek(8,1)
f.close()

# find the firing points:
Vind = np.where(V==0)[0]
k = 0
print Vind
plt.figure(figsize=(16,16), dpi = 80)
plt.subplot(2,1,1)
plt.plot(np.arange(len(V))*0.5,V,label = 'V')
k = 0
while (k < len(Vind)):
    plt.axvline(Vind[k]*0.5, color='red')
    k += 4
plt.subplot(2,1,2)
plt.plot(np.arange(len(I))*0.5,I,label = 'I')
k = 0
while (k < len(Vind)):
    plt.axvline(Vind[k]*0.5, color='red')
    k += 4
#plt.show()
plt.savefig('V_I.png')
plt.close()
