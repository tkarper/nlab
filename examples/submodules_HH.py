import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

# Discretization parameters
nStell = 100		# Number of stellate cells
dt = 0.01

# Connection strengths
I = 0.0			# External exitatory input
th2s = 1.0		# Theta oscillation input to stellates
s2in = 1.0		# Stellate to interneuron
in2s = -1.0		# Interneuron to stellate
hd2s = 3.0		# HD cell to stellate
cvar.Neuron_Ack_gh = 0


def get_neuron_item(neurArr, varname):
	arr = np.zeros(neurArr.shape)
	for i  in range(neurArr.size):
		arr[i] = type(neurArr[i]).__swig_getmethods__.get(varname, None)(neurArr[i])
	return arr



print('Initializing...')



# Create stellate cells
#Stellates = np.array([Neuron_HH() for _ in range(0,nStell)])
#cvar.Neuron_HH_alpha = 2
#cvar.Neuron_HH_beta = 2
Stellates = np.array([Neuron_Ack() for _ in range(0,nStell)])
for i in range(0,nStell):
	Stellates[i].VP += 2*(random.random()-1)
	Stellates[i].I = I


# Create theta oscillator and link to stellates
Theta = np.array([Neuron_Osc(20, 3, 1)])	# Period, duration, strength
connect_one_to_many(Theta[0], Stellates, th2s)
#for i in range(0,nStell):
#	Stellates[i].I = I

		
# Create interneurons and link to one another
nIntNeuro = int(0.05*nStell)		# Number of interneurons
IntNeuros = np.array([Neuron_Ack() for _ in range(0,nIntNeuro)])
connect_many_to_many(IntNeuros, IntNeuros, in2s)

# Link interneurons to stellate cells
nSPerIN = int(0.6*nStell)	# Number of stellates per submodule
## Uniformly random sample
#for i in range(0,nIntNeuro):
#	submodInd = random.sample(xrange(nStell), nSPerIN)
#	submod = Stellates[submodInd]
#	connect_many_to_one(submod, IntNeuros[i], s2in)
#	connect_one_to_many(IntNeuros[i], submod, in2s)
## Lateral normal distribution
distr = np.random.normal(0.0, 0.2, (nStell, nIntNeuro))
distr = distr
for i in range(0,nIntNeuro):
	# Convert normally distributed numbers to indices in [0, nStell)
	submodInd = np.round(i*nStell/float(nIntNeuro) + nStell*distr[...,i]).astype(int)
	# Compute indices modulo nStell, avoiding negative index values
	submodInd = np.mod(nStell + np.mod(submodInd, nStell), nStell)
#	print(submodInd)
	# Extract the submodule and create connections to and from interneuron i
	submod = Stellates[submodInd]
	connect_many_to_one(submod, IntNeuros[i], s2in)
	connect_one_to_many(IntNeuros[i], submod, in2s)

# Create head cells
nHead = 2		# Number of head cells
#Heads = np.array([Neuron() for _ in range(0,nHead)])
Heads = np.array([Neuron_Osc(10, 3, 1), Neuron_Osc(5, 0, 0)])

# Link HD cells to stellates
Nh2sConns = 1	# Number of head to stellate connections per stellate cell
for i in range(0,nStell):
	headInd = random.sample(xrange(nHead), Nh2sConns)
	connect_many_to_one(Heads[headInd], Stellates[i], hd2s)
#	for j in range(0,Nh2sConns):
#		headInd = random.randrange(nHead)
#		Stellates[i].connect(Heads[headInd], hd2s)


print('Initialization finished.')


t = 0.0
m = 0
plotInd = 0
plt.figure()
plt.ion()


tHist = []
hist = []

# MAIN TIMELOOP
while(True):
	t = t+dt
	m = m+1
	
#	# Update HD cells
#	direction = 0.0
##	direction = fmod(2*pi + pi/8.0*sin(t/10.0), 2*pi)		# Direction in radians, from x-axis
#	for i in range(0, nHead):
#		b = direction/(2*pi)
#		c = 0.1
#		Heads[i].s = max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
	
	# Update neural network
	stepNetwork(Theta, t, dt)
	stepNetwork(Stellates, t, dt)
	stepNetwork(IntNeuros, t, dt)
	stepNetwork(Heads, t, dt)
	
	# Check if stellates have fired
	for i in range(0, 20):
		if Stellates[i].V>0 and Stellates[i].V < Stellates[i].VP:
				tHist.append(t)
				hist.append(i)
	
	# Move to next timestep
	updateNetwork(Stellates)
	updateNetwork(IntNeuros)
	updateNetwork(Theta)
	updateNetwork(Heads)
	
	# Turn off external input for an initial phase
	if (t < 10):
		Theta[0].sp = 0
		for i in range(0,nHead):
			Heads[i].sp = 0
	
	# Plot data
	plotInterval = 10.0
	if(fmod(t, plotInterval) < dt):
		print('t=%f'%t)
#		print('Theta:sp = %f'%Theta[0].sp)
		
		plt.clf()
		sp = get_neuron_item(Stellates, 'V')
		plt.subplot(2,3,1)
		plt.plot(sp, 'gs')
		plt.ylim((-100,100))
		plt.title('Stellate memb. pot.')
		
		plt.subplot(2,3,2)
		sp = get_neuron_item(Stellates, 'msyn')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Stellate syn. pot.')
		
		plt.subplot(2,3,3)
		sp = get_neuron_item(IntNeuros, 'V')
		plt.plot(sp, 'gs')
		plt.ylim((-100,100))
		plt.title('Interneuron memb. pot.')
		
		plt.subplot(2,3,4)
#		sp = get_neuron_item(Heads, 's')
#		plt.plot(sp, 'bo')
#		plt.ylim((0,1))
#		plt.title('Head cells')
		plt.plot(tHist, hist, 'rs')
		plt.xlim((t-400,t))
		plt.title('Firing history')
		
		
		plt.subplot(2,3,6)
		sp = get_neuron_item(IntNeuros, 's')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Interneuron syn. pot.')
		
		
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
