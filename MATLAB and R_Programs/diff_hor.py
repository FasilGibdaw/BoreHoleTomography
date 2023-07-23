import numpy as np

def diff_hor(row, col):
    # Function to compute 1st order difference along the horizontal
    L = np.diag(np.ones(col)) - np.diag(np.ones(col - 1), -1)
    L1 = np.kron(np.eye(row), L)
    return L1
