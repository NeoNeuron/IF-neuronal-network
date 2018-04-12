#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

dt = 0.05

def g(t, t0, s, tau):
    g_val = np.zeros(len(t))
    t0_ind = int(t0/dt)
    g_val[t0_ind::] = s*np.exp(-(t[t0_ind::]-i)/tau)
    return g_val

t = np.arange(0,80,dt)
#gl = np.arange(0.01, 0.15, 0.005)
#gl = np.array([0.005])
gl = 0.005
ee = 14.0/3.0
# ge = np.zeros(len(t))
# ti = np.arange(0,40,1000)
# s = np.arange(0.001, 0.051, 0.001)
s = 0.005
# tau = np.arange(1, 3.1, 0.1)
tau = 2
v_ini = np.arange(0.9,1,0.01)
I_max = np.zeros(len(v_ini))
I_min = np.zeros(len(v_ini))
counter = 0
ge = s*np.exp(-t/tau)
for j in v_ini:
	v = np.ones(len(t))*j
	for i in range(len(t) - 1):
		v[i + 1] = v[i] + dt * (-gl*v[i] - ge[i]*(v[i] - ee))
		if v[i + 1] >= 1:
			v[i + 1] = 0
	I = -gl * v - ge*(v - ee)
	I_min[counter] = I.min()
	I_max[counter] = I.max()
	counter += 1
        plt.figure(figsize = (16,6))
        plt.subplot(1,2,1)
        plt.plot(t, v, label = 'V')
        plt.subplot(1,2,2)
        plt.plot(t, I, label = 'I')
	plt.legend()
	plt.show()

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
# plt.plot(gl, v_tau)
plt.plot(v_ini, I_max, label = 'I_max')
plt.plot(v_ini, I_min, label = 'I_min')
plt.xlabel('Initial Potential Level')
plt.ylabel('Total synaptic current')
plt.legend()
plt.show()
