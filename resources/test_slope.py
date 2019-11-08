import numpy as np
import matplotlib.pyplot as plt
from eksternlab import height, slope, curvature
from scipy.constants import g


x = [  0.0,   0.2,   0.4,   0.6,   0.8,   1.0,   1.2,   1.4]
y = [0.472, 0.374, 0.313, 0.280, 0.284, 0.314, 0.373, 0.450]

h = height(x, y)

x_interp = np.linspace(0.0, 1.0, 50)

alpha = slope(h, x_interp)
kappa = curvature(h, x_interp)
plt.plot(h, x_interp)
plt.ylabel('$slope$'), plt.xlabel('$x$')
plt.show()