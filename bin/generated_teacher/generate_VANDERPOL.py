import matplotlib.pyplot as plt
import numpy as np
nt = 60000
dt = 0.03
mu = 1
def fx(x, y):
    return y
def fy(x, y):
    return mu * (1 - x * x) * y - x
dimteach = 2
x = np.zeros(nt)
y = np.zeros(nt)
x[0] = 0.001
y[0] = 0.01
for i in range(nt - 1):
    x[i + 1] = x[i] + dt * fx(x[i], y[i])
    y[i + 1] = y[i] + dt * fy(x[i], y[i])
plt.plot(x, y)
plt.show()
step = len(x)
file = open('VANDERPOLdt=0_03', 'w')
file.write(str(len(x)) + ' ' + str(2))
file.write('\n')
for i in range(len(x)):
    file.write(str(x[i]) + ' ')
    file.write(str(y[i]) + '\n')
file.close()
