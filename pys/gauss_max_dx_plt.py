#!/usr/bin/python
# This script aims to study the dependence of maximum mutual information of the parematric model of Gaussian random variables;
# The results can be ploted as functions of mutual interacting strength, for both theoretical value and numerical results.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
from scipy.optimize import curve_fit
def func_fit(x, a):
    return a*x**2
a = 0.5
b = 0.5
axy = 1
target_ind = 0
bins = np.arange(0.1,0.51,0.01)
# import data
mis = np.loadtxt(sys.argv[1])
# rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)
rho2 = axy**2/(1+a**2+axy**2)
theory = -0.5*np.log(1-rho2)
dev = abs(mis - theory)
# fitting
popt, pcov = curve_fit(func_fit, bins, dev)
print popt
# R2
var2 = sum((dev - dev.mean())**2)
R = 1 - sum((dev - func_fit(bins, popt[0]))**2)/var2

fig = plt.figure()
plt.plot(bins, dev, '.', label = 'data') 
plt.plot(bins, func_fit(bins, popt[0]), label = 'fitting')
plt.title('R = %.6f' % R)
ax = plt.gca()
#ax.xaxis.get_major_formatter().set_powerlimits((0,1))
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.xlabel('Binsize')
plt.ylabel('Deviation(I_true - I_obs)')
plt.legend()
plt.grid()
if len(sys.argv) == 2:
    plt.show()
elif len(sys.argv) == 3:
    plt.savefig(sys.argv[2])
    plt.close()
