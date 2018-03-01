#!/usr/bin/python
# Peakcock test, two-dimensional KS test for two emperical data;
import numpy as np
import matplotlib.pyplot as plt
import sys
# # test data;
# D1 = np.empty((2,10000))
# D2 = np.zeros((2,10000))
# D1[0] = np.random.normal(0, 1, 10000)
# D1[1] = np.random.normal(0, 1, 10000)
# D2[0] = np.random.normal(0, 2, 10000)
# D2[1] = np.random.normal(0, 2, 10000)
# # print D1
# # print D2

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
                cdf[i][0] = H[i][0] * A
            else:
                cdf[i][j] = cdf[i][j - 1] + H[i][j] * A
        if i != 0:
            cdf[i] += cdf[i - 1]
    return cdf

def Peacock2d(X, Y, bin_num_1, bin_num_2):
    # Find the range of histograms;
    max1 = X[0].max() if X[0].max() >= Y[0].max() else Y[0].max()
    min1 = X[0].min() if X[0].min() <= Y[0].min() else Y[0].min()
    max2 = X[1].max() if X[1].max() >= Y[1].max() else Y[1].max()
    min2 = X[1].min() if X[1].min() <= Y[1].min() else Y[1].min()
    # Calculate histograms;
    H1, edges1, edges2 = np.histogram2d(X[0], X[1], bins = [bin_num_1, bin_num_2], range = [[min1, max1], [min2, max2]], normed = True)
    H2, edges1, edges2 = np.histogram2d(Y[0], Y[1], bins = [bin_num_1, bin_num_2], range = [[min1, max1], [min2, max2]], normed = True)
    C1ff = cdf2d(H1, edges1, edges2, orient = (False, False))
    C1ft = cdf2d(H1, edges1, edges2, orient = (False, True))
    C1tf = cdf2d(H1, edges1, edges2, orient = (True, False))
    C1tt = cdf2d(H1, edges1, edges2, orient = (True, True))

    C2ff = cdf2d(H2, edges1, edges2, orient = (False, False))
    C2ft = cdf2d(H2, edges1, edges2, orient = (False, True))
    C2tf = cdf2d(H2, edges1, edges2, orient = (True, False))
    C2tt = cdf2d(H2, edges1, edges2, orient = (True, True))

    diff_ff = abs(C1ff - C2ff).max()
    diff_tf = abs(C1tf - C2tf).max()
    diff_ft = abs(C1ft - C2ft).max()
    diff_tt = abs(C1tt - C2tt).max()

    diff_max = max(diff_ff, diff_ft, diff_tf, diff_tt)
    p = np.exp(-D1.shape[1]*diff_max**2)
    return diff_max, p

# Settings for histograms;
bin_num_1 = 2
bin_num_2 = 100

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
D, p = Peacock2d(D1, D2, bin_num_1, bin_num_2)
print('D = %.6f'%D)
print('p-value = %.6f'%p)
