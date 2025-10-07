import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style



plt.figure(1)
ax = plt.axes(projection = '3d')
x1 = np.arange(-5, 5, 0.25)
y1 = np.arange(-5, 5, 0.25)
X1, Y1 = np.meshgrid(x1, y1)
Z1 = np.sin(X1) * np.cosh(Y1)
ax.plot_surface(X1, Y1, Z1, cmap = 'Spectral')
ax.set_title('3D Surface Plot')

plt.figure(2)
ax = plt.axes(projection = '3d')
x2 = np.arange(-5, 5, 0.25)
y2 = np.arange(-5, 5, 0.25)
X2, Y2 = np.meshgrid(x1, y1)
Z2 = np.sin(X2) * np.sin(Y2)
ax.plot_surface(X2, Y2, Z2, cmap = 'Spectral')
ax.set_title('3D Surface Plot')
plt.show()