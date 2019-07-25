import matplotlib.pyplot as plt
import numpy as np
nt = 60000
dt = 0.03
mu = 1
def teacherchaos(x, y, J, a, d, eps, beta):
    for i in range(len(x) - 1):
        y[i + 1] = y[i] + eps * (x[i] - J)
        x[i + 1] = x[i] + x[i] * (x[i] - a) * (1 - x[i]) - beta * (x[i] > d) - y[i]
dimteach = 2

x = np.zeros(nt)
y = np.zeros(nt)
x[0] = .1
y[0] = .2
J1 = 0.327; eps1 = 0.00793; beta1 = 0.19547; d1 = 0.5012; a1 = 0.25
teacherchaos(x, y, J1, a1, d1, eps1, beta1)
#x *= 2
#y *= 10
#%matplotlib
#plt.scatter(x, y, s = .1)
#plt.show()
plt.scatter(x, y, s = .1)
plt.show()
step = len(x)
file = open('Neuron_Chaos=J1 = 0.327; eps1 = 0.00793; beta1 = 0.19547; d1 = 0.5012; a1 = 0.25', 'w')
file.write(str(len(x)) + ' ' + str(2))
file.write('\n')
for i in range(len(x)):
    file.write(str(x[i]) + ' ')
    file.write(str(y[i]) + '\n')
file.close()
