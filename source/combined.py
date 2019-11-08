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

def Heun(fv, fx, fn, v_0, x_0, h, N): # fA = function, A_0 = initial value, h = step size, N = no. of points
    # Init arrays for simulation variables
    t = np.zeros(N)
    v = np.zeros(N)
    x = np.zeros(N)
    norm = np.zeros(N)
    alpha = np.zeros(N)
    kappa = np.zeros(N)

     # Set IVs
    v[0] = v_0
    x[0] = x_0

    # Calculation
    for n in range(N-1):
        # increment t by h
        t[n+1] = t[n] + h

        # update alpha and kappa at current x
        alpha[n] = exlab.slope(lane_h, x[n])
        kappa[n] = exlab.curvature(lane_h, x[n])

        # calculate next x
        x[n+1] = x[n] + h*v[n]

        # calculate next v
        v_star = v[n] + h*fv(v[n], t[n], alpha[n])
        v[n+1] = v[n] + 0.5*h*( fv(v[n], t[n], alpha[n]) + fv(v_star, t[n], alpha[n]) )

        # calculate next norm
        norm_star = norm[n] + h*fn(norm[n], t[n], alpha[n], 1/kappa[n])
        norm[n+1] = norm[n] + 0.5*h*( fn(norm[n], t[n], alpha[n], 1/kappa[n]) + fn(norm_star, t[n], alpha[n], 1/kappa[n]) )

    return v, x, norm, alpha, t

def plot_heun(fv, fx, fn, v_0, x_0, h, N):
    # Set figure params and plot
    graph_params = {'figure.figsize': (16, 12), 'axes.grid': True,
                'lines.linewidth': 3, 'font.size': 14}
    plt.rcParams.update(graph_params)

    # Run simulation
    v_Heun, x_Heun, n_Heun, alpha, t = Heun(fv, fx, fn, v_0, x_0, h, N)

    # find abs(v(t))
    for i in range(len(v_Heun)):
        v_Heun[i] = abs(v_Heun[i])
    
    # plot v(t)
    plt.subplot(2, 1, 1)
    plt.title("RK2 numerical solution of v(t)")
    plt.plot(t, v_Heun)
    plt.ylabel('$v(t)$'), plt.xlabel('$t$')

    # plot x(t)
    plt.subplot(2, 1, 2)
    plt.title("RK2 numerical solution of x(t)")
    plt.plot(t, x_Heun)
    plt.ylabel('$x(t)$'), plt.xlabel('$t$')

    plt.show()

def plot_exp():
    # data is of the form  t   x   y   v   a   r
    lane_data = np.loadtxt("/home/christopher/Documents/3. semester/TFY4115/resources/data.txt")
    data_t = lane_data[:, 0]
    data_x = lane_data[:, 1]
    data_y = lane_data[:, 2]
    data_v = lane_data[:, 3]
    data_a = lane_data[:, 4]
    
    # v(t)
    plt.subplot(2, 1, 1)
    plt.title("Experimental data for v(t)")
    plt.plot(data_t, data_v)
    plt.ylabel('$v(t)$')

    # x(t)
    plt.subplot(2, 1, 2)
    plt.title("Experimental data for x(t)")
    plt.plot(data_t, data_x)
    plt.ylabel('$x(t)$'), plt.xlabel('$t$')
    plt.show()

def plot_comp(fv, fx, fn, v_0, x_0, h, N):
    # Set figure params and plot
    graph_params = {'figure.figsize': (16, 12), 'axes.grid': True,
                'lines.linewidth': 3, 'font.size': 14}
    plt.rcParams.update(graph_params)

    # Run simulation
    v_Heun, x_Heun, n_Heun, alpha, t = Heun(fv, fx, fn, v_0, x_0, h, N)

    # find abs(v(t))
    for i in range(len(v_Heun)):
        v_Heun[i] = abs(v_Heun[i])

    # data is of the form  t   x   y   v   a   r
    lane_data = np.loadtxt("/home/christopher/Documents/3. semester/TFY4115/resources/data.txt")
    data_t = lane_data[:, 0]
    data_x = lane_data[:, 1]
    data_y = lane_data[:, 2]
    data_v = lane_data[:, 3]
    data_a = lane_data[:, 4]
    
    # v(t)
    plt.subplot(2, 1, 1)
    plt.title("Experimental data for v(t)")
    plt.plot(data_t, data_v)
    plt.ylabel('$v(t)$')

    plt.subplot(2, 1, 2)
    plt.title("RK2 numerical solution of v(t)")
    plt.plot(t, v_Heun)
    plt.ylabel('$v(t)$'), plt.xlabel('$t$')
    plt.show()

    # x(t)
    plt.subplot(2, 1, 1)
    plt.title("Experimental data for x(t)")
    plt.plot(data_t, data_x)
    plt.ylabel('$x(t)$')
   
    plt.subplot(2, 1, 2)
    plt.title("RK2 numerical solution of x(t)")
    plt.plot(t, x_Heun)
    plt.ylabel('$x(t)$'), plt.xlabel('$t$')
    plt.show()


if __name__ == '__main__':
    h = 0.001
    N = 30000
    x_0 = 0
    v_0 = 0
    
    # Function for dv(t)/dt
    fv = lambda v, t, alpha : (7*g*np.sin(alpha) - 7*k*v/m) / 5
    # Function for dx(t)/dt
    fx = lambda x, t, alpha : 0
    # Function for N(t)
    fn = lambda n, v, alpha, r : - m*g*np.cos(alpha) + m*v**2 / r

    plot_comp(fv, fx, fn, v_0, x_0, h, N)
