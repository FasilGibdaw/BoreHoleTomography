import numpy as np
import matplotlib.pyplot as plt

dx = 1  # step along the horizontal
dy = 1  # step along the vertical
x = np.arange(0, 11, dx)
y = np.arange(0, 11, dy)
nx, ny = len(x) - 1, len(y) - 1  # grid nx by ny
xmin, xmax = x.min(), x.max()  # end points of the horizontal
ymin, ymax = y.min(), y.max()  # end points of the vertical
rangex, rangey = (xmax - xmin), (ymax - ymin)

side1 = 5  # Number of receivers along side 1
Rxy1 = np.column_stack(((xmin + (np.arange(1, side1 + 1) * rangex) / (side1 + 1)), np.full(side1, ymax)))
side2 = 5  # Number of receivers along side 2
Rxy2 = np.column_stack((np.full(side2, xmax), ymin + (np.arange(1, side2 + 1) * rangey) / (side2 + 1)))
Rxy = np.vstack((Rxy1, Rxy2))  # Receiver positions along side 1 and 2

side3 = 5  # Number of transmitters along side 3
Txy1 = np.column_stack((np.full(side3, xmin), ymin + (np.arange(1, side3 + 1) * rangey) / (side3 + 1)))
side4 = 0  # No transmitters along side 4
Txy2 = np.column_stack(((xmin + (np.arange(1, side4 + 1) * rangex) / (side4 + 1)), np.full(side4, ymin)))
Txy = np.vstack((Txy1, Txy2))

b = 0  # initializing number of rays
N_ray = (side1 + side2) * (side3 + side4)  # number of rays from transmitter-receiver pair
ray = np.zeros((N_ray, 4))  # initializing ray matrix
for t in range(side3 + side4):
    for r in range(side1 + side2):
        ray[b, 0] = Txy[t, 0]  # x components of the transmitters position
        ray[b, 1] = Txy[t, 1]  # y components of the transmitters position
        ray[b, 2] = Rxy[r, 0]  # x components of the receivers position
        ray[b, 3] = Rxy[r, 1]  # y components of the receivers position
        b += 1

# Graphics
plt.figure(figsize=(8, 6))
# Grid plot
for j in range(len(y)):
    plt.plot(x, y[j] * np.ones(len(x)), 'k', linewidth=1)
for k in range(len(x)):
    plt.plot(x[k] * np.ones(len(y)), y, 'k', linewidth=1)

# Plot position of the transmitters and receivers
plt.plot(Rxy[:, 0], Rxy[:, 1], 'o', markersize=8, markerfacecolor='r', label='Receivers')
plt.plot(Txy[:, 0], Txy[:, 1], 'rs', markersize=8, markerfacecolor='g', label='Transmitters')

# Plot ray path between transmitters and receivers
for i in range(N_ray):
    plt.plot([ray[i, 0], ray[i, 2]], [ray[i, 1], ray[i, 3]],'k')

plt.xlabel('Position (arbitrary units)')
plt.ylabel('Position (arbitrary units)')
plt.legend()
plt.savefig('fig1.png')
plt.show()
