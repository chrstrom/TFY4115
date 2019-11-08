
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import eksternlab as exlab
from scipy.constants import g

# Main simulation loop happens here
def sim(h, N, v_0, x_0):
    v = np.zeros(N); v[0] = v_0
  
    x = np.zeros(N); x[0] = x_0
    t = np.zeros(N)
    norm = np.zeros(N)
    fric = np.zeros(N)
    alpha = np.zeros(N)

    # Eulers shitty method
    for n in range(N-1):
        t[n+1] = t[n] + h
        alpha[n] = exlab.slope(lane_h, x[n])

        x[n+1] = x[n] + h*v[n]*np.cos(alpha[n])
        v[n+1] = v[n] + 0.5*h*fv(v[n], alpha[n])
        
        norm[n] = fnorm(v[n], x[n], alpha[n])
        fric[n] = ffric(v[n], x[n], alpha[n])

    return t, v, x, norm, fric

# Data functions; retrieval and fixing
def get_exp_data(filename, FORMAT = False):
    
    if(not FORMAT):
        # data is of the form  t   x   y   v   a   r
        exp_data = np.loadtxt(filename)
        exp_data_t = exp_data[:, 0]
        exp_data_x = exp_data[:, 1]
        exp_data_y = exp_data[:, 2]
        exp_data_v = exp_data[:, 3]
        exp_data_a = exp_data[:, 4]
    else:
        # data is of the form  t   x   y   r   v   a
        exp_data = np.loadtxt(filename)
        exp_data_t = exp_data[:, 0]
        exp_data_x = exp_data[:, 1]
        exp_data_y = exp_data[:, 2]
        exp_data_v = exp_data[:, 4]
        exp_data_a = exp_data[:, 5] 

    return exp_data_t, exp_data_v, exp_data_x

def fix_exp_data(exp_data_t, exp_data_v, exp_data_x, FORMAT = False):
    # DIfferent parameters if format specified True
    if(not FORMAT):
        exp_x_offset = -0.05
        exp_v_offset = 0.1
    else:
        exp_x_offset = -0.05
        exp_v_offset = 0
        
    # Because v_exp is relative to the horizon and not the track, we need to change this
    for n in range(len(exp_data_v)):
        alpha = exlab.slope(lane_h, exp_data_x[n])
        exp_data_v[n] *= np.cos(alpha)
        
    # Offset data to start at the origin
    exp_t_0 = exp_data_t[0]
    exp_x_0 = exp_data_x[0]
    exp_v_0 = exp_data_v[0]
    exp_data_v[0] = 0

    for i in range(len(exp_data_t)):
        exp_data_t[i] -= exp_t_0
        exp_data_x[i] -= exp_x_offset
        exp_data_v[i] -= exp_v_0 - exp_v_offset

    return exp_data_v, exp_data_x
def fix_sim_data(v, norm, fric):
    # Set v = |v| for all n, to match exp data
    for i in range(len(v)):
        v[i] = abs(v[i])

    # Remove last data point for norm
    norm[-1] = norm[-2]

    # Remove last data point for fric
    fric[-1] = fric[-2]

    return v, norm, fric

# Plotting functions
def set_plot_params():
      # Set figure params and plot
    graph_params = {'figure.figsize': (16, 12), 'axes.grid': True,
                'lines.linewidth': 3, 'font.size': 14}
    plt.rcParams.update(graph_params)

def plot_comp(t1, t2, f1, f2, labelf1, labelf2, xlabel, ylabel, title = ""):    
    f1_plt, = plt.plot(t1, f1)
    f2_plt, = plt.plot(t2, f2)
    plt.legend((f1_plt, f2_plt), (labelf1, labelf2))
    plt.title(title)
    plt.xlabel(xlabel), plt.ylabel(ylabel)
    plt.show()

if __name__ == '__main__':

    # Constants for simulation
    k = 0.00505067561487144             # Calculated with k_calc.py
    m = 0.0298                          # Mass of the ball           

    h = 0.001                           # Time step                             
    N = 35000                           # Amount of steps
    v_0 = 0                             # Initial velocity
    x_0 = 0                             # Initial position

    data_file = "/home/christopher/Documents/3. semester/TFY4115/resources/data.txt"
    shifted_data_file = "/home/christopher/Documents/3. semester/TFY4115/resources/shifted_axis_data.txt"
# Functions
    fv = lambda v, alpha : ( 7 * g * np.sin(alpha) - 7 * k * v*np.cos(alpha) / m ) / 5
    fnorm = lambda v, x, alpha : m * g * np.cos(alpha) - np.cos(alpha + np.pi) * m * v**2 * exlab.curvature(lane_h, x)
    ffric = lambda v, x, alpha: m * g * np.sin(alpha) - m * ( 7 * g * np.sin(alpha) - 7 * k * v / m ) / 5

# Lane data, measurements from lab
    lane_data_x = [0.000, 0.200, 0.400, 0.600, 0.800,  1.00,  1.20,  1.40]
    lane_data_y = [0.472, 0.374, 0.313, 0.280, 0.284, 0.314, 0.373, 0.450]

    lane_h = exlab.height(lane_data_x, lane_data_y)

# Run simulation and fix sim data
    t_sim, v_sim, x_sim, norm_sim, fric_sim = sim(h, N, v_0, x_0)           # All functions are in the scope of sim
    v, norm, fric = fix_sim_data(v_sim, norm_sim, fric_sim)
   
# Get experimental data from file and fix said data
    t_exp, v_exp, x_exp = get_exp_data(data_file)
    t_exp2, v_exp2, x_exp2 = get_exp_data(shifted_data_file, True)

    v_exp, x_exp = fix_exp_data(t_exp, v_exp, x_exp)
    v_exp2, x_exp2 = fix_exp_data(t_exp2, v_exp2, x_exp2, True)

# Plots
    # Update plot parameters
    set_plot_params()

    # Experimental and simulation data for v(t)
    plot_comp(t_sim, t_exp, v_sim, v_exp, "Simulering", "Eksperimentell data", "t [s]", "v [m/s]")

    # Experimental and simulation data for x(t)
    plot_comp(t_sim, t_exp, x_sim, x_exp, "Simulering", "Eksperimentell data", "t [s]", "x [m]")

    # Normal force on the ball from the track and friction between the ball and the track
    plot_comp(t_sim, t_sim, fric_sim, norm_sim, "Friksjon", "Normalkraft", "t [s]", "Kraft [N]")

    # Experimenting with different values for k
    # This runs the simulation an extra two times; consider commenting out unless needed
        # k = 0.05
        # t_1, v_1, x_1, norm_1, fric_1 = sim(h, N, v_0, x_0) 

        # k = 0.0005
        # t_2, v_2, x_2, norm_2, fric_2 = sim(h, N, v_0, x_0) 

        # plot_comp(t_1, t_2, x_1, x_2, "k = 0.05", "k = 0.0005", "t [s]", "x [m]")

    # Plot comparison between normal data and shifted data:
    plot_comp(t_exp, t_exp2, v_exp, v_exp2, "Normal data", "Shifted data", "t [s]", "v [m/s]")
