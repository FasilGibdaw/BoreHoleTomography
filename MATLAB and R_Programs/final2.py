import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from diff_hor import diff_hor
from diff_ver import diff_ver
import cv2
# Load the input image and preprocess it
#A = plt.imread('endi5.png').astype(float)
A = cv2.imread('endi5.png',cv2.IMREAD_UNCHANGED)
ima = -(A[:, :, 0] - 256)  # true model
dx, dy = 1, 1
ima = ima[::dx, ::dy]  # sampled model
ny, nx = ima.shape

# Grid points along the x and y axis
x = np.arange(nx + 1)
y = np.arange(ny + 1)

# End points of the x and y axis
xmin, xmax = x.min(), x.max()
ymin, ymax = y.min(), y.max()

# Define receiver and transmitter positions
side1 = 5  # Number of receivers along side 1 at the surface of the ground
rangex = xmax - xmin
rangey = ymax - ymin
Rxy1 = np.array([(xmin + ((np.arange(1, side1 + 1) * rangex) / (side1 + 1)), ymax * np.ones(side1))]).T
side2 = 5  # Number of receivers along side 2 (borehole 2)
Rxy2 = np.array([(xmax * np.ones(side2), ymin + ((np.arange(1, side2 + 1) * rangey) / (side2 + 1)))], dtype=float).T
Rxy = np.vstack((Rxy1, Rxy2))  # Receiver positions along side 1 and 2

side3 = 5  # Number of transmitters along side 3 (borehole 1)
Txy1 = np.array([(xmin * np.ones(side3), ymin + ((np.arange(1, side3 + 1) * rangey) / (side3 + 1)))], dtype=float).T
side4 = 0
Txy2 = np.array([(xmin + ((np.arange(1, side4 + 1) * rangex) / (side4 + 1)), ymin * np.ones(side4))]).T
Txy = np.vstack((Txy1, Txy2))

# Construct rays with corresponding x and y coordinates
b = 0
N_ray = (side1 + side2) * (side3 + side4)
ray = np.zeros((N_ray, 4))  # 4 columns for p1(x1, y1) and p2(x2, y2) of the ray
for t in range(side3 + side4):
    for r in range(side1 + side2):
        ray[b, 0] = Txy[t, 0]  # x components of p1
        ray[b, 1] = Txy[t, 1]  # y components of p1
        ray[b, 2] = Rxy[r, 0]  # x components of p2
        ray[b, 3] = Rxy[r, 1]  # y components of p2
        b += 1

pixel = np.arange(nx * ny) + 1  # number of pixels in one vector

# Calculating the theory matrix for the forward problem to simulate the measurements
Nr = 3000  # number of discretization along each ray
G = np.zeros((N_ray, len(pixel)))  # initializing the theory matrix
for n in range(N_ray):
    d = np.sqrt((ray[n, 0] - ray[n, 2]) ** 2 + (ray[n, 1] - ray[n, 3]) ** 2)  # length of a ray
    dr = d / Nr  # step length of a ray
    drx = (ray[n, 2] - ray[n, 0]) / Nr  # x coordinate of step length
    dry = (ray[n, 3] - ray[n, 1]) / Nr  # y coordinate of step length
    xp = ray[n, 0] + np.arange(1, Nr + 1) * drx  # number of x coordinates along a ray
    yp = ray[n, 1] + np.arange(1, Nr + 1) * dry  # number of y coordinates along a ray
    ind_x = np.floor(((xp - xmin) / rangex) * (nx - 1)).astype(int)
    ind_y = np.floor(((yp - ymin) / rangey) * (ny - 1)).astype(int)
    index = (ind_y - 1) * nx + ind_x
    dd, cc = np.unique(index, return_counts=True)
    l_ind = cc * dr
    G[n, dd] = G[n, dd] + l_ind

# Generating the simulated measurement m with error E
E = 13.3 * np.random.randn(N_ray)  # assumed error in the measurement
m = G @ ima.flatten() + E  # the measurement

# Discretization for the inverse problem
dx1, dy1 = 0.7, 0.7  # discretization step different from dx and dy from above
x1 = np.arange(0, nx, dx1)  # grid points along the x axis
y1 = np.arange(0, ny, dy1)  # grid points along the y axis
nx1, ny1 = len(x1) - 1, len(y1) - 1
x1min, x1max = x1.min(), x1.max()  # end points of the x axis
y1min, y1max = y1.min(), y1.max()  # end points of the y axis
rangex1, rangey1 = x1max - x1min, y1max - y1min
pixel_1 = np.arange(nx1 * ny1) + 1  # number of pixels in one vector

# Calculating the theory matrix for the inversion problem
G1 = np.zeros((N_ray, len(pixel_1)))  # initializing the theory matrix
for n in range(N_ray):
    d = np.sqrt((ray[n, 0] - ray[n, 2]) ** 2 + (ray[n, 1] - ray[n, 3]) ** 2)  # length of a ray
    dr = d / Nr  # step length of a ray
    drx = (ray[n, 2] - ray[n, 0]) / Nr  # x coordinate of step length
    dry = (ray[n, 3] - ray[n, 1]) / Nr  # y coordinate of step length
    xp = ray[n, 0] + np.arange(1, Nr + 1) * drx  # number of x coordinates along a ray
    yp = ray[n, 1] + np.arange(1, Nr + 1) * dry  # number of y coordinates along a ray
    ind_x = np.floor(((xp - x1min) / rangex1) * (nx1 - 1)).astype(int)
    ind_y = np.floor(((yp - y1min) / rangey1) * (ny1 - 1)).astype(int)
    index = (ind_y - 1) * nx1 + ind_x
    dd, cc = np.unique(index, return_counts=True)
    l_ind = cc * dr
    G1[n, dd] = G1[n, dd] + l_ind

# Constructing the prior (difference prior) and MAP estimate
e = 0.2 * E
sigma = np.var(e)  # assumed measurement error covariance
L1=diff_ver(nx1,ny1);#function to compute difference matrix along horizontal(L1)  
L2=diff_hor(nx1,ny1);#function to compute difference matrix along vertical(L2)
LL = np.vstack((L2, L1))  # difference matrix (horizontal + vertical)
C = 1.2 * np.random.randn(2 * nx1 * ny1)  # virtual measurement error (uncertainty in the difference)
cc = np.var(C) * np.eye(2 * nx1 * ny1)  # covariance of the uncertainty in the difference
# Xm = np.linalg.inv((G1.T @ np.linalg.inv(sigma * np.eye(len(m))) @ G1) + LL.T @ np.linalg.inv(cc) @ LL) @ (G1.T @ np.linalg.inv(sigma * np.eye(len(m))))
# tolerance = 1e-8
# U, s, VT = np.linalg.svd((G1.T @ np.linalg.inv(sigma * np.eye(len(m))) @ G1) + LL.T @ np.linalg.inv(cc) @ LL)
# S_inv = np.zeros((len(s), G1.T.shape[0]))
# S_inv[:len(s), :len(s)] = np.diag(1.0 / s)
# Xm = VT.T @ S_inv @ U.T @ (G1.T @ np.linalg.inv(sigma * np.eye(len(m))) @ m)
# Xm = Xm.reshape(ny1, nx1)  # Reshape Xm to 2D
Xm = np.linalg.inv((G1.T @ np.linalg.inv(sigma * np.eye(len(m))) @ G1) + LL.T @ np.linalg.inv(cc) @ LL) @ (G1.T @ np.linalg.inv(sigma * np.eye(len(m))))
XmapN = Xm @ m  # MAP estimate
XmapN = XmapN.reshape(ny1, nx1)  # back to the original dimension (2D)
ima = ima.reshape(ny, nx)  # back to the original dimension (2D) of the assumed true model

# Graphics
plt.figure(figsize=(15, 5))
plt.subplot(131)
plt.pcolormesh(ima)
plt.title('(a) True model (625 pixels)')

plt.subplot(132)
plt.pcolormesh(ima)
for i in range(0, N_ray):  # DN is the step size to select rays to be plotted
    plt.plot([ray[i, 0], ray[i, 2]], [ray[i, 1], ray[i, 3]], 'c')
plt.axis('equal')
plt.title('(b) True model with rays')

plt.subplot(133)
plt.pcolormesh(XmapN)
plt.title('(c) MAP estimate (1225 pixels)')

plt.tight_layout()
plt.show()
