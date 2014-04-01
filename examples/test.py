import sys
sys.path.append('../')

from nlab import *
from math import *
import numpy as np
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

nInter = 2
nSPIN = 5	# Stellates Per InterNeuron
nStell = nSPIN*nInter
stell = np.array([Neuron_Ack() for _ in range(0,nStell)])
inter = np.array([Neuron_Ack() for _ in range(0,nInter)])

for i in range(0,nStell):
	stell[i].I = I #* (1 + (i/nSPIN)*1.0) #I*(1+0.1*random.random())
#	stell[i].connect(theta, th2s*(1+i*0.1/nStell))
#	stell[i].VP += random.random()*10
for i in range(nInter):
	submod = stell[i*nInter:(i+1)*nInter]
	connect_many_to_one(submod, inter[i], s2in)
	connect_one_to_many(inter[i], submod, in2s)

connect_one_to_many(theta, stell, th2s)

t = 0.0
m = 0


tHist = []
hist = [[] for _ in range(nStell)]
intHist = [[] for _ in range(nInter)]
thetaHist = []

# MAIN TIMELOOP
while(t<500):
	t = t+dt
	m = m+1
	
	# Update neural network
	stepNetwork(inter, t, dt)
	theta.step(t,dt,0)
	stepNetwork(stell, t, dt)
	updateNetwork(inter)
	theta.update()
	updateNetwork(stell)
	
	if (m%10 == 0):
		tHist.append(t)
		for i in range(nStell):
			hist[i].append(stell[i].V)
		for i in range(nInter):
			intHist[i].append(inter[i].V)
		thetaHist.append(theta.s)


nPlot = nStell+2
plt.figure()
plt.ion()
for i in range(nStell):
	plt.subplot(nPlot,1,i+1)
	plt.plot(tHist,hist[i])
	plt.xlim((0, tHist[-1]))
for i in range(nInter):
	plt.subplot(nPlot,1,nStell+1+i)
	plt.plot(tHist,intHist[i],'r')
	plt.xlim((0, tHist[-1]))
plt.subplot(nPlot,1,nStell+nInter)
plt.plot(tHist,thetaHist,'g')
plt.xlim((0, tHist[-1]))
