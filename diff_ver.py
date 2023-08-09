import numpy as np

def diff_ver(row, col):
    # Function to compute 1st order difference along the vertical
    c = np.zeros(row * col)
    c[0] = 1
    c[col] = -1
    L2 = np.zeros((row*col, row*col))
    for i in range(row*col):
        L2[i, i:] = c[:row*col - i]
    return L2
