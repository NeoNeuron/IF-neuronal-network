#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

def g(t, t0):
    gg = np.zeros(len(t))
    dt = 0.1
    tau = 2
    s = 0.005
    t0_ind = int(t0/dt)
    gg[t0_ind::] = s*np.exp(-(t[t0_ind::]-i)/tau)
    return gg

dt = 0.1
t = np.arange(0,1000,dt)
v = np.zeros(len(t))
ge = np.zeros(len(t))
ti = np.arange(0,1000,100)

for i in ti:
    ge += g(t, i)

for i in range(9999):
    v[i + 1] = v[i] + dt * (-0.05*v[i] - ge[i]*(v[i] - 14.0/3.0))
    if v[i + 1] > 1:
        v[i + 1] = 0
I = -0.05 * v - ge*(v - 14.0/3.0)

# print len(t)
# print len(ge)
plt.plot(t, v)
# plt.show()
# plt.hist(ge,200)
plt.show()
