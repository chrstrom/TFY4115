# Implementation of Heun's method (RK2)
import numpy as np
import matplotlib.pyplot as plt
import eksternlab as exlab
from scipy.constants import g

# CONSTANTS
k = 0.005067267482262274
m = 0.0298

# Measure in 8 points
lane_data_x = [  0.0,   0.2,   0.4,   0.6,   0.8,   1.0,   1.2,   1.4]
lane_data_y = [0.472, 0.374, 0.313, 0.280, 0.284, 0.314, 0.373, 0.450]

lane_h = exlab.height(lane_data_x, lane_data_y)

def Heun(fv, v_0, h, N): # fA = function, A_0 = initial value, h = step size, N = no. of points
    # Init arrays for simulation variables
    t = np.zeros(N)
    v = np.zeros(N)

    N_MAX = 10000

    # test plot for alpha over one interval
    x_interp = np.linspace(0.0, 1.0, N_MAX+1)
    alpha = exlab.slope(lane_h, x_interp)

     # Set IVs
    v[0] = v_0

    counter = 0
    # Calculation
    for n in range(N-1):
        # increment t by h
        t[n+1] = t[n] + h

        if (counter%2 != 0):
            # calculate next v
            v_star = v[n] + h*fv(v[n], t[n], alpha[n%N_MAX])
            v[n+1] = v[n] + 0.5*h*( fv(v[n], t[n], alpha[n%N_MAX]) + fv(v_star, t[n], alpha[n%N_MAX]) )
        else:
            v_star = v[n] + h*fv(v[n], t[n], alpha[N_MAX - n%N_MAX])
            v[n+1] = v[n] + 0.5*h*( fv(v[n], t[n], alpha[N_MAX - n%N_MAX]) + fv(v_star, t[n], alpha[N_MAX - n%N_MAX]) )

        if (n % 10000 == 0):
            counter += 1

    return v, t

def plotHeun(fv, v_0, h, N):
    # Set figure params and plot
    graph_params = {'figure.figsize': (16, 6), 'axes.grid': True,
                'lines.linewidth': 3, 'font.size': 14}
    plt.rcParams.update(graph_params)

    # Run simulation
    v_Heun, t = Heun(fv, v_0, h, N)
    # plot v(t)
    plt.title("RK2 numerical solution of v(t)")
    plt.plot(t, v_Heun)
    plt.ylabel('$v(t)$'), plt.xlabel('$t$')
    plt.show()


if __name__ == '__main__':
    h = 0.0005
    N = 60000
    x_0 = 0
    v_0 = 0
    
    # Function for dv(t)/dt
    fv = lambda v, t, alpha : (7*g*np.sin(alpha) - 7*k*v/m) / 5

    # Plot solution v(t)
    plotHeun(fv, v_0, h, N)
