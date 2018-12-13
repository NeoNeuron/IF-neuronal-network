#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import struct as st
import subprocess

f = open("data/spike/singleSpike.bin", 'rb')
shape = st.unpack('QQ', f.read(16))
dat0 = np.array(st.unpack('b'*shape[1], f.read(1*shape[1])))
f.close()
dat0 = dat0[:-1]

bin_size = 0.005

#f = open("data/tmp/singleI.bin", 'rb')
f = open("data/lfp/singleI_con.bin", 'rb')
shape = st.unpack('QQ', f.read(16))
dat1 = np.array(st.unpack('d'*shape[1], f.read(8*shape[1])))
f.close()
dat1 = dat1[1:]
#bins1 = np.arange(min(dat1), max(dat1) + bin_size, bin_size)
dat1_1 = dat1[dat0==1] 
#print(len(dat1_1))

f = open("data/lfp/singleI_nocon.bin", 'rb')
shape = st.unpack('QQ', f.read(16))
dat2 = np.array(st.unpack('d'*shape[1], f.read(8*shape[1])))
f.close()
dat2=dat2[1:]
if min(dat1) < min(dat2):
    MIN = min(dat1)
else:
    MIN = min(dat2)
if max(dat1) > max(dat2):
    MAX = max(dat1)
else:
    MAX = max(dat2)

bins = np.arange(MIN, MAX + bin_size, bin_size)
#print(bins)

fig, axes = plt.subplots(2, 3, sharex = True, figsize = (20,12), dpi = 80)
########################################
h1 = axes[0,0].hist(dat1, bins = bins, density = True)
#axes[0,0].set_xlabel('Y( Local Field Potential )')
axes[0,0].set_ylabel('P(Y)')
axes[0,0].set_title('Connected neuron pair')
axes[0,0].grid(linestyle='--')
########################################
h2 = axes[0,1].hist(dat2, bins = bins, density = True)
#axes[0,1].set_xlabel('Y( Local Field Potential )')
axes[0,1].set_ylabel('P(Y)')
axes[0,1].set_title('Independent neuron pair')
axes[0,1].grid(linestyle='--')
########################################
axes[0,2].plot(h1[1][:-1], h1[0]-h2[0])
#axes[0,2].set_xlabel("Y")
axes[0,2].set_ylabel("Probability difference")
axes[0,2].grid(linestyle='--')
########################################
h4 = axes[1,0].hist(dat1_1, bins = bins, density = True)
axes[1,0].set_xlabel('Y( Local Field Potential )')
axes[1,0].set_ylabel('P(Y)')
axes[1,0].set_title("P(y|x=1)")
axes[1,0].grid(linestyle='--')
########################################
axes[1,1].plot(h1[1][:-1], h4[0]-h1[0])
axes[1,1].set_xlabel('Y( Local Field Potential )')
axes[1,1].set_ylabel("Probability difference")
axes[1,1].set_title("P(y|x=1)-P(y))")
axes[1,1].grid(linestyle='--')
########################################
axes[1,2].plot(h1[1][2:], -1*np.diff(h1[0], n=1))
axes[1,2].set_xlabel('Y( Local Field Potential )')
axes[1,2].set_ylabel("Probability difference")
axes[1,2].set_title("P(y-dy)-P(y))")
axes[1,2].grid(linestyle='--')
########################################

plt.savefig('lfp_hist')
plt.close()
#plt.show()
