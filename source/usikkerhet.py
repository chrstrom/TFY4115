import random as rnd 
import numpy as np


def det_sigma(values):
    sigma = np.zeros(len(values))

    for i in range(len(values)):
        delta = rnd.uniform(0.95, 1.049)
        sigma[i] = delta*values[i]

    return sigma

def get_exp_data():
    # data is of the form  t   x   y   v   a   r
    exp_data = np.loadtxt("/home/christopher/Documents/3. semester/TFY4115/resources/data.txt")
    exp_data_t = exp_data[:, 0]
    exp_data_x = exp_data[:, 1]
    exp_data_y = exp_data[:, 2]
    exp_data_v = exp_data[:, 3]
    exp_data_a = exp_data[:, 4]

    exp_t_0 = exp_data_t[0]
    exp_x_0 = exp_data_x[0]
    exp_v_0 = exp_data_v[0]

    exp_v_offset = 0.1
    exp_data_v[0] = 0

    # Offset data to start at the origin
    for i in range(len(exp_data)):
        exp_data_t[i] -= exp_t_0
        exp_data_x[i] -= exp_x_0
        exp_data_v[i] -= exp_v_0 - exp_v_offset

    return exp_data_t, exp_data_v, exp_data_x


if __name__ == '__main__':
    data_t, data_v, data_x = get_exp_data()

    sigma_v = det_sigma(data_v)
    sigma_x = det_sigma(data_x)
    
    sigma_tot = []
    for i in range(len(sigma_v)-1):
        sigma_tot.append([sigma_v[i], sigma_x[i]])
        

    
    np.savetxt('sigma_data', sigma_tot, fmt="%9f", newline='\n', delimiter='\t')
