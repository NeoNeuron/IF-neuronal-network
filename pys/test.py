#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

dt = 0.05

def g(t, t0, s, tau):
    g_val = np.zeros(len(t))
    t0_ind = int(t0/dt)
    g_val[t0_ind::] = s*np.exp(-(t[t0_ind::]-i)/tau)
    return g_val

t = np.arange(0,160,dt)
gl = np.arange(0.01, 0.15, 0.005)
ee = 14.0/3.0
# ge = np.zeros(len(t))
# ti = np.arange(0,40,1000)
# s = np.arange(0.001, 0.051, 0.001)
s = 0.005
# tau = np.arange(1, 3.1, 0.1)
tau = 2
v_max = np.zeros(len(gl))
I_min = np.zeros(len(gl))
v_tau = np.zeros(len(gl))
for j in range(len(gl)):
    ge = s*np.exp(-t/tau)
    v = np.zeros(len(t))
    for i in range(len(t) - 1):
        v[i + 1] = v[i] + dt * (-gl[j]*v[i] - ge[i]*(v[i] - ee))
        if v[i + 1] >= 1:
            v[i + 1] = 0
    v_max[j] = v.max()
    v_max_ind = v.argmax()
    v_diff = abs(v - v.max()*np.exp(-1))
    v_diff_min_ind = v_diff.argmin()
    v_tau[j] = (v_diff_min_ind - v_max_ind) * dt
    I = -gl[j] * v - ge*(v - ee)
    I_min[j] = I.min()
    # plt.plot(t, v, label = str(gl[j]))


# plt.plot(s, v_max, label = 'data')
# v_fit = np.polyfit(s, v_max, 1)
# print v_fit
# plt.plot(s, v_fit[0]*s + v_fit[1], label = 'fit')
# plt.xlabel('Strength')
# plt.ylabel('IPSP(mV)')
# plt.title('IPSP = %.2f*s + %.2e'%(v_fit[0],v_fit[1]))
# plt.legend()
# plt.show()
# plt.hist(ge,200)
plt.plot(gl, v_tau)
plt.legend()
plt.show()
