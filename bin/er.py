import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('error')

#plt.plot(data[:, 1], data[:, 0])
#plt.show()

#plt.plot(data)
#plt.show()
data2 = np.loadtxt('generated_teacher/van')
out = data2[0:50000, :] - data

plt.plot(data2[:, 0], data2[:, 1])
plt.plot(out[40000:50000, 0], out[40000:50000, 1])

plt.show()