#!/usr/bin/python
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys

def func_fit1(x, a):
    return a*x
def func_fit2(x, a):
    return a*x**2
dx = np.arange(0.001, 0.0101, 0.0002)
# import data
arr = np.loadtxt('mi_max_dx_1e8.csv')
dev = abs(arr-arr[0])
# linear fitting
popt1, pcov1 = curve_fit(func_fit1, dx[1:], dev[1:])
print popt1
# square fitting
popt2, pcov2 = curve_fit(func_fit2, dx[1:], dev[1:])
print popt2
# R2
var2 = sum((dev[1:] - dev[1:].mean())**2)
R1 = 1 - sum((dev[1:] - func_fit1(dx[1:], popt1[0]))**2)/var2
R2 = 1 - sum((dev[1:] - func_fit2(dx[1:], popt2[0]))**2)/var2
fig = plt.figure()
plt.plot(dx[1:], dev[1:], '.', label = 'data')
#plt.plot(dx[1:], func_fit1(dx[1:], popt1[0]), label = 'linear fitting')
plt.plot(dx[1:], func_fit2(dx[1:], popt2[0]), label = 'square fitting')
#plt.title('R1 = %.4f, R2 = %.4f' %(R1, R2))
plt.title('R = %.4f' %R2)
ax = plt.gca()
ax.xaxis.get_major_formatter().set_powerlimits((0,1))
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.xlabel('Binsize')
plt.ylabel('Deviation(I_true - I_obs)')
plt.grid()
plt.legend()
if len(sys.argv) == 1:
    plt.show()
elif len(sys.argv) == 2:
    plt.savefig(sys.argv[1])
    plt.close()
