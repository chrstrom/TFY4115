# Implementation of Heun's method (RK2)

import numpy as np
import matplotlib.pyplot as plt
from eksternlab import height, slope, curvature
from scipy.constants import g
# CONSTANTS
k = 0.0005
m = 0.004

# Set figure params and plot
graph_params = {'figure.figsize': (16, 6), 'axes.grid': True,
             'lines.linewidth': 2, 'font.size': 14}
plt.rcParams.update(graph_params)

### We can use this to map the experimental data for x, y, t ###
# Import lane data for t, x, and y
# lane_data = np.loadtxt("C:/Users/Christopher/Desktop/TFY4115/bane.txt")
# t_data = lane_data[:, 0]
# lane_data_x = lane_data[:, 1]
# lane_data_y = lane_data[:, 2]

# Measure in 8 points
lane_data_x = []
lane_data_y = []
lane_h = height(lane_data_x, lane_data_y)

def Heun(f, v_0, h, N): # f = function, x_0 = initial value, h = step size, N = no. of points

    # Init arrays for t and v
    t = np.zeros(N)
    v = np.zeros(N)

     # Set IVs
    v[0] = v_0

    # Splitting the entire slope into N discrete elements
    x_interp = np.linspace(0.0, 1.0, N)
    # y_interp = lane_h(x_interp)

    # Calculation
    for n in range(N-1):
        # Update theta
        theta = slope(lane_h, x_interp)
      
        v_star = v[n] + h*f(v[n], t[n], theta)
        v[n+1] = v[n] + 0.5*h*(f(v[n], t[n], theta) + f(v_star, t[n], theta))
        
    return v, t

def plotHeun(f, v_0, h, N, y_label, x_label):
    # Run simulation
    v_Heun, t = Heun(f, v_0, h, N)

    # Plot results
    plt.plot(t,v_Heun)
    plt.ylabel(y_label), plt.xlabel(x_label)
    plt.show()

if __name__ == '__main__':
    h = 0.05
    N = 2000
    v_0 = 0
    
    # Function for dv(t)/dt
    def dvdt(v, t, theta):
        return (5*g*np.sin(theta) - 5*k*v/m) / 3
    
    # Function for N(t)
    def norm(v, t, theta):
        return theta * t

    # Plot solution v(t)
    plotHeun(dvdt, v_0, h, N, '$v(t)$', '$t$')
    plotHeun(norm, v_0, h, N, '$v(t)$', '$t$')