import numpy as np
import matplotlib.pyplot as plt
#sin

def sinus(dt, T0, T1):
    pi = np.pi
    Time = np.arange(T0, T1, dt)
    data = np.sin(10 * np.pi * Time)
    return data

T0 = 0
T1 = 50
dt = 0.001
data = sinus(dt, T0, T1)
file = open('sin(10pi)_dt_0_001', 'w')
file.write(str(len(data)) + ' ' + str(1))
file.write('\n')
for i in range(len(data)):
    file.write(str(data[i]) + '\n')
file.close()