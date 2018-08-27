#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit

dat = np.loadtxt(sys.argv[1], delimiter=',')

plt.plot(dat[:,0], dat[:,1], '.')
plt.xlabel('Connecting probability')
plt.ylabel('Maximum mutual information')
# plot settings:
xmajorLocator   = MultipleLocator(1) 
#xmajorFormatter = FormatStrFormatter('$10^{%d}$') # set formats of x axis
xminorLocator   = MultipleLocator(0.1) 

ymajorLocator   = MultipleLocator(1e-4) 
#ymajorFormatter = FormatStrFormatter('%1.1f') # set formats of y axis
yminorLocator   = MultipleLocator(0.2e-4) 

ax = plt.gca()
#ax.xaxis.set_major_formatter(xmajorFormatter)
#ax.xaxis.set_major_locator(xmajorLocator)
#ax.xaxis.set_minor_locator(xminorLocator)
#
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
#ax.yaxis.set_major_locator(ymajorLocator)
#ax.yaxis.set_minor_locator(yminorLocator)

# square fit:
def func(x, a, b):
    return a*x**2 + b

popt, pcov = curve_fit(func, dat[:,0], dat[:,1])
#print popt
#print pcov
plt.plot(np.arange(0.1,0.7,0.05), func(np.arange(0.1,0.7,0.05), popt[0], popt[1]), label = 'fitting curve')
plt.legend()
plt.grid()
plt.savefig('lfp_mi-p.png')
#plt.show()
