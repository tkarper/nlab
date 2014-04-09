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
I = 8.0			# External exitatory input
th2s = 0 #7.0
in2s = 5
in2in = 0#5
s2in = 10
nInter = 1
nSPerIN = 10	# Stellates Per InterNeuron



def get_neuron_item(neurArr, varname):
	arr = np.zeros(neurArr.shape)
	for i  in range(neurArr.size):
		arr[i] = type(neurArr[i]).__swig_getmethods__.get(varname, None)(neurArr[i])
	return arr



nStell = nSPerIN*nInter
stell = np.array([Neuron_Stel() for _ in range(0,nStell)])
inter = np.array([Neuron_Pyr() for _ in range(0,nInter)])


# Create submodule indices 'submod' and set initial data
submod = [[] for _ in range(nInter)]
for i in range(nInter):
	submod[i] = stell[i*nSPerIN:(i+1)*nSPerIN]
#	submod[i] = stell[random.sample(xrange(nStell), nSPerIN)]
	for j in range(nSPerIN):
		submod[i][j].V += random.random()*20
		submod[i][j].I = I #* (1 + i*5)
#		stell[i].connect(theta, th2s*(1+i*0.1/nStell))

# Connect stellates <-> interneurons <-> interneurons
for i in range(nInter):
	connect_many_to_one(submod[i], inter[i], s2in)
	connect_one_to_many(inter[i], submod[i], in2s)
connect_many_to_many(inter, inter, in2in)




t = 0.0
m = 0


tHist = []
sVHist = [[] for _ in range(nStell)]
sSHist = [[] for _ in range(nStell)]
intVHist = [[] for _ in range(nInter)]




# MAIN TIMELOOP
updateNetwork(inter)
updateNetwork(stell)
while(t<200):
	t = t+dt
	m = m+1
	
	# Update neural network
	stepNetwork(inter, t, dt)
	stepNetwork(stell, t, dt)
	updateNetwork(inter)
	updateNetwork(stell)
	
	
	if (m%10 == 0):
		tHist.append(t)
		for i in range(nStell):
			sVHist[i].append(stell[i].V)
			sSHist[i].append(stell[i].s)
		for i in range(nInter):
			intVHist[i].append(inter[i].V)


# PLOT DATA
plotStelS = True
plt.figure()
plt.ion()
if plotStelS:
	nPlot = 2*nStell+nInter+1
else:
	nPlot = nStell+nInter+1
ctr = 1
for i in range(nInter):
	for j in range(nSPerIN):
		k = i*nSPerIN+j
		plt.subplot(nPlot,1,ctr)
		ctr += 1
		plt.plot(tHist,sVHist[k], 'r')
		plt.xlim((0, tHist[-1]))
		if plotStelS:
			plt.subplot(nPlot,1,ctr)
			ctr += 1
			plt.plot(tHist,sSHist[k], 'g')
			plt.xlim((0, tHist[-1]))
	plt.subplot(nPlot,1,ctr)
	ctr += 1
	plt.plot(tHist,intVHist[i],'b')
	plt.xlim((0, tHist[-1]))
#plt.subplot(nPlot,1,nStell+nInter+1)
#plt.plot(tHist,thetaHist,'g')
#plt.xlim((0, tHist[-1]))
