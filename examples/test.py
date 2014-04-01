import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

# Discretization parameters
dt = 0.01

# Connection strengths
I = 2.0			# External exitatory input
cvar.Neuron_Ack_gh = 0


def get_neuron_item(neurArr, varname):
	arr = np.zeros(neurArr.shape)
	for i  in range(neurArr.size):
		arr[i] = type(neurArr[i]).__swig_getmethods__.get(varname, None)(neurArr[i])
	return arr



nStell = 5
inter = Neuron_Ack()
stell = np.array([Neuron_Ack() for _ in range(0,nStell)])
for i in range(0,nStell):
	stell[i].I = I

connect_many_to_one(stell, inter, 2)
connect_one_to_many(inter, stell, -2.0)

t = 0.0
m = 0


tHist = []
hist1 = []
hist2 = []

# MAIN TIMELOOP
while(t<100):
	t = t+dt
	m = m+1
	
	# Update neural network
	inter.step(t,dt,0)
	stepNetwork(stell, t, dt)
	inter.update()
	updateNetwork(stell)
	
	if (m%10 == 0):
		tHist.append(t)
		hist1.append(stell[0].V)
		hist2.append(inter.V)


plt.figure()
plt.ion()
plt.subplot(1,2,1)
plt.plot(tHist,hist1)
plt.subplot(1,2,2)
plt.plot(tHist,hist2)
