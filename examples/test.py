import numpy as np
import matplotlib.pyplot as plt


def mollifyDiracs(x, t, sigma):
	x = np.transpose(np.atleast_2d(x))
	y = np.exp(-(x-t)**2 / sigma**2)
	return y.sum(0)

x = np.array([1, 2, 4, 5, 8])
t = np.linspace(0,10)
y = mollifyDiracs(x, t, 1)

print(t)
print(y)

plt.figure()
plt.ion()
plt.plot(t, y)
