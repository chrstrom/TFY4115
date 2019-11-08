import numpy as np
from scipy.constants import g

m = 0.0298

# Read data from files, put in lane_data
# data is of the form  t   x   y   v   a   r
lane_data = np.loadtxt("/home/christopher/Documents/3. semester/TFY4115/resources/data.txt")
data_v = lane_data[:, 3]

for i in range(len(data_v)):
    data_v[i] = data_v[i]**2            # squaring all the data

y_start = lane_data[0, 2]               # start y-val
y_end = lane_data[-1, 2]                # end y-val

dx = lane_data[1, 0] - lane_data[0, 0]  # Determined by tracker dataset
I = np.trapz([data_v], dx = dx)
k = m*g*(y_start - y_end) / I

print("I = " + str(I))
print("k = " + str(k))