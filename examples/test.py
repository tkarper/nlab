import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

cvar.Neuron_Ack_gh = 0

# Discretization parameters
dt = 0.01

# Connection strengths
I = 0.0			# External exitatory input
th2s = 7.0
in2s = -2
s2in = 2


def get_neuron_item(neurArr, varname):
	arr = np.zeros(neurArr.shape)
	for i  in range(neurArr.size):
		arr[i] = type(neurArr[i]).__swig_getmethods__.get(varname, None)(neurArr[i])
	return arr


theta = Neuron_Osc(20, 3, 1)

nStell = 5
nInter = 2
stell = np.array([Neuron_Ack() for _ in range(0,nStell)])
inter = Neuron_Ack()
for i in range(0,nStell):
	stell[i].I = I #*(1+i*0.1/nStell) #I*(1+0.1*random.random())
#	stell[i].connect(theta, th2s*(1+i*0.1/nStell))
#	stell[i].VP += random.random()*10

connect_many_to_one(stell, inter, s2in)
connect_one_to_many(inter, stell, in2s)
connect_one_to_many(theta, stell, th2s)

t = 0.0
m = 0


tHist = []
hist = [[] for _ in range(nStell)]
intHist = []
thetaHist = []

# MAIN TIMELOOP
while(t<500):
	t = t+dt
	m = m+1
	
	# Update neural network
	inter.step(t,dt,0)
	theta.step(t,dt,0)
	stepNetwork(stell, t, dt)
	inter.update()
	theta.update()
	updateNetwork(stell)
	
	if (m%10 == 0):
		tHist.append(t)
		for i in range(nStell):
			hist[i].append(stell[i].V)
		intHist.append(inter.V)
		thetaHist.append(theta.s)


nPlot = nStell+2
plt.figure()
plt.ion()
for i in range(nStell):
	plt.subplot(nPlot,1,i+1)
	plt.plot(tHist,hist[i])
	plt.xlim((0, tHist[-1]))
plt.subplot(nPlot,1,nStell+1)
plt.plot(tHist,intHist,'r')
plt.xlim((0, tHist[-1]))
plt.subplot(nPlot,1,nStell+2)
plt.plot(tHist,thetaHist,'g')
plt.xlim((0, tHist[-1]))
