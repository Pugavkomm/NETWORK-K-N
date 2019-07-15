import numpy as np
import matplotlib.pyplot as plt

def periodic_impuls(T0, T1, dt):
    Time = np.arange(T0, T1, dt)
    data = np.sin(10 * np.pi * Time)
    t = -1/5/4 + 1/5
    f_in = np.zeros(len(Time))
    for i in range(len(Time)):
        if t >= 1/5:
            f_in[i] = 1;
            t = 0.
        t += dt
        
    plt.plot(data)
    plt.plot(f_in)
    plt.show()
    file = open('f_in_for_sin_10_dt_0_001', 'w')
    file.write(str(1) + '\n')
    for i in range(len(Time)):
        file.write(str(f_in[i]) + '\n')
    file.close()

T0 = 0
T1 = 50
dt = 0.001
periodic_impuls(T0, T1, dt)