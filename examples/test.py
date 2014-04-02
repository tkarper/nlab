import sys
sys.path.append('../')

from nlab import *
from math import *
import numpy as np
import random
import matplotlib.pyplot as plt



# Discretization parameters
dt = 0.01

# Connection strengths
I = 1.5			# External exitatory input
th2s = 0 #7.0
in2s = 0 #-2
in2in = 0 #-10.0
s2in = 3


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
inter = np.array([Neuron_Ack() for _ in range(0,nInter)])\


for i in range(nStell):
	stell[i].I = I

submod = [[] for _ in range(nInter)]
for i in range(nInter):
	submod[i] = stell[i*nSPIN:(i+1)*nSPIN]
	connect_many_to_one(submod[i], inter[i], s2in)
	connect_one_to_many(inter[i], submod[i], in2s)

for i in range(nInter):
	connect_one_to_many(theta, submod[i], th2s*(1+i*0.1))
#	for j in range(nSPIN):
#		submod[i][j].I = I * (1 + i*0.3)
#		stell[i].connect(theta, th2s*(1+i*0.1/nStell))
#		submod[i][j].VP += random.random()*20
connect_many_to_many(inter, inter, in2in)
#connect_one_to_many(theta, stell, th2s)

t = 0.0
m = 0


tHist = []
sHist = [[] for _ in range(nStell)]
intHist = [[] for _ in range(nInter)]
thetaHist = []



# MAIN TIMELOOP
while(t<600):
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
			sHist[i].append(stell[i].V)
		for i in range(nInter):
			intHist[i].append(inter[i].V)
		thetaHist.append(theta.s)


nPlot = nStell+nInter+1
plt.figure()
plt.ion()
for i in range(nInter):
	for j in range(nSPIN):
		k = i*nSPIN+j
		plt.subplot(nPlot,1,k+i+1)
		plt.plot(tHist,sHist[k], 'r')
		plt.xlim((0, tHist[-1]))
	plt.subplot(nPlot,1,(i+1)*nSPIN+i+1)
	plt.plot(tHist,intHist[i],'b')
	plt.xlim((0, tHist[-1]))
plt.subplot(nPlot,1,nStell+nInter+1)
plt.plot(tHist,thetaHist,'g')
plt.xlim((0, tHist[-1]))
