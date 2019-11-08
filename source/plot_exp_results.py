
import numpy as np
import matplotlib.pyplot as plt

# Set figure params and plot
graph_params = {'figure.figsize': (16, 6), 'axes.grid': True,
            'lines.linewidth': 1, 'font.size': 14}
plt.rcParams.update(graph_params)

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
plt.ylabel('$v(t)$'), plt.xlabel('$t$')

# a(t)
plt.subplot(2, 1, 2)
plt.title("Experimental data for a(t)")
plt.plot(data_t, data_a)
plt.ylabel('$a(t)$'), plt.xlabel('$t$')
plt.show()