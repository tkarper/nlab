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
I = 5.0			# External exitatory input
tMax = 1000



def get_neuron_item(neurArr, varname):
	arr = np.zeros(neurArr.shape)
	for i  in range(neurArr.size):
		arr[i] = type(neurArr[i]).__swig_getmethods__.get(varname, None)(neurArr[i])
	return arr



neur = np.array([Neuron_Stel(), Neuron_Inter()])
nNeur = neur.size



neur[1].EK = -90
neur[1].ENa = 70
neur[1].gL = 0.5
neur[1].gK = 11
neur[1].gNa = 52
neur[1].gNap = 0.5


t = 0.0
m = 0


fireHist = [[] for _ in range(nNeur)]
tFireHist = [[] for _ in range(nNeur)]
tHist = []
sVHist = [[] for _ in range(nNeur)]
sSHist = [[] for _ in range(nNeur)]




# MAIN TIMELOOP
#updateNetwork(neur)
while(t<tMax):
	t = t+dt
	m = m+1
	
	newI = I*t/tMax
	neur[0].I = newI
	neur[1].I = newI
	
	# Update neural network
	stepNetwork(neur, t, dt)
	updateNetwork(neur)
	
	
	# Check if stellates have fired
	for i in range(0, nNeur):
		if neur[i].isFiring:
			fireHist[i].append(neur[i].V)	
			tFireHist[i].append(t)	
	
	if (m%10 == 0):
		tHist.append(t)
		for i in range(nNeur):
			sVHist[i].append(neur[i].V)
			sSHist[i].append(neur[i].s)


# Compute frequencies
freq = [[] for _ in range(nNeur)]
for i in range(nNeur):
	nFire = max(len(fireHist[i])-1, 0)
	freq[i] = np.zeros(nFire)
	for j in range(nFire):
		freq[i][j] = 1000/(tFireHist[i][j+1] - tFireHist[i][j])


# PLOT DATA
plotStelS = True
plt.figure()
plt.ion()
if plotStelS:
	nPlot = 3*nNeur
else:
	nPlot = 2*nNeur
ctr = 1
for i in range(nNeur):
	plt.subplot(nPlot,1,ctr);	ctr += 1
	plt.plot(tHist,sVHist[i], 'r')
	plt.xlim((0, tHist[-1]))
	if plotStelS:
		plt.subplot(nPlot,1,ctr);		ctr += 1
		plt.plot(tHist,sSHist[i], 'g')
		plt.xlim((0, tHist[-1]))
	plt.subplot(nPlot,1,ctr);		ctr += 1
	plt.plot(tFireHist[i][1:], freq[i], 'bo-')
	plt.xlim((0, tHist[-1]))
	plt.ylabel('Freq')

