#!/usr/bin/python
# Peakcock test, two-dimensional KS test for two emperical data;
import numpy as np
import matplotlib.pyplot as plt
import sys

# Settings for histograms;
x_bin_num = 2
y_bin_num = 100

# # test data;
# D1 = np.empty((2,10000))
# D2 = np.zeros((2,10000))
# D1[0] = np.random.normal(0, 1, 10000)
# D1[1] = np.random.normal(0, 1, 10000)
# D2[0] = np.random.normal(0, 2, 10000)
# D2[1] = np.random.normal(0, 2, 10000)
# # print D1
# # print D2

# import original data; Suppose D1 and D2 are 2 by N arrays;
boolmat = np.genfromtxt(sys.argv[1], delimiter = ',')
doublemat = np.genfromtxt(sys.argv[2], delimiter = ',')
bool_ind_1 = int(sys.argv[3])
bool_ind_2 = int(sys.argv[4])
double_ind_1 = int(sys.argv[5])
double_ind_2 = int(sys.argv[6])

D1 = np.concatenate(([boolmat[bool_ind_1]], [doublemat[double_ind_1]]))
D2 = np.concatenate(([boolmat[bool_ind_2]], [doublemat[double_ind_2]]))
D1 = np.delete(D1, -1, 1)
D2 = np.delete(D2, -1, 1)
# print D1.shape
# print D2
# Find the range of histograms;
xmax = D1[0].max() if D1[0].max() >= D2[0].max() else D2[0].max()
xmin = D1[0].min() if D1[0].min() <= D2[0].min() else D2[0].min()
ymax = D1[1].max() if D1[1].max() >= D2[1].max() else D2[1].max()
ymin = D1[1].min() if D1[1].min() <= D2[1].min() else D2[1].min()

# Calculate histograms;
H1, xedges, yedges = np.histogram2d(D1[0], D1[1], bins = [x_bin_num, y_bin_num], range = [[xmin, xmax], [ymin, ymax]], normed = True)
H2, xedges, yedges = np.histogram2d(D2[0], D2[1], bins = [x_bin_num, y_bin_num], range = [[xmin, xmax], [ymin, ymax]], normed = True)
def cdf2d(H, xedges, yedges, orient = (False, False)):
    # Prerequirest: H is normed;
    cdf = np.empty(H.shape)
    xstep = xedges[1] - xedges[0]
    ystep = yedges[1] - yedges[0]
    A = xstep * ystep
    # orientation selection:
    if orient[0]:
        H = np.flipud(H)
    if orient[1]:
        H = np.fliplr(H)
    # Start Calculation;
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            if j == 0:
                # cdf[i][0] = H[i][0] * A
            else:
                cdf[i][j] = cdf[i][j - 1] + H[i][j] * A
        if i != 0:
            cdf[i] += cdf[i - 1]
    return cdf

C1ff = cdf2d(H1, xedges, yedges, orient = (False, False))
C1ft = cdf2d(H1, xedges, yedges, orient = (False, True))
C1tf = cdf2d(H1, xedges, yedges, orient = (True, False))
C1tt = cdf2d(H1, xedges, yedges, orient = (True, True))

C2ff = cdf2d(H2, xedges, yedges, orient = (False, False))
C2ft = cdf2d(H2, xedges, yedges, orient = (False, True))
C2tf = cdf2d(H2, xedges, yedges, orient = (True, False))
C2tt = cdf2d(H2, xedges, yedges, orient = (True, True))

diff_ff = abs(C1ff - C2ff).max()
diff_tf = abs(C1tf - C2tf).max()
diff_ft = abs(C1ft - C2ft).max()
diff_tt = abs(C1tt - C2tt).max()

diff_max = max(diff_ff, diff_ft, diff_tf, diff_tt)
print('D = %.6f'%diff_max)
print('p-value = %.6f'%(np.exp(-D1.shape[1]*diff_max**2)))
