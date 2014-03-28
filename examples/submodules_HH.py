import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

# Discretization parameters
nStell = 50		# Number of stellate cells
dt = 0.01

# Connection strengths
I = 0.5			# External exitatory input
s2in = 1.0		# Stellate to interneuron
in2s = -1.0		# Interneuron to stellate
hd2s = 2.0		# HD cell to stellate



print('Initializing...')



# Create stellate cells
Stellates = np.array([Neuron_HH() for _ in range(0,nStell)])
cvar.Neuron_HH_alpha = 2
cvar.Neuron_HH_beta = 2

# Create theta oscillator and link to stellates
Theta = np.array([Neuron_Osc(10, 1, 1)])	# Period, duration, strength
connect_one_to_many(Theta[0], Stellates, I)
#for i in range(0,nStell):
#	Stellates[i].I = I
		
# Create interneurons and link to one another
nIntNeuro = 5	# Number of interneurons
IntNeuros = np.array([Neuron_HH() for _ in range(0,nIntNeuro)])
connect_many_to_many(IntNeuros, IntNeuros, in2s)

# Link interneurons to stellate cells
nSPerIN = 30	# Number of stellates per submodule
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
	print(submodInd)
	# Extract the submodule and create connections to and from interneuron i
	submod = Stellates[submodInd]
	connect_many_to_one(submod, IntNeuros[i], s2in)
	connect_one_to_many(IntNeuros[i], submod, in2s)

# Create head cells
nHead = 2		# Number of head cells
Heads = np.array([Neuron() for _ in range(0,nHead)])
#Heads = np.array([Neuron_Osc(5, 1, 1), Neuron_Osc(5, 0, 0)])

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


# MAIN TIMELOOP
while(True):
	t = t+dt
	m = m+1
	
	# Update HD cells
	direction = 0.0
#	direction = fmod(2*pi + pi/8.0*sin(t/10.0), 2*pi)		# Direction in radians, from x-axis
	for i in range(0, nHead):
		b = direction/(2*pi)
		c = 0.1
		Heads[i].s = max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
	
	# Update neural network
	stepNetwork(Theta, t, dt)
	stepNetwork(Stellates, t, dt)
	stepNetwork(IntNeuros, t, dt)
	stepNetwork(Heads, t, dt)
	updateNetwork(Stellates)
	updateNetwork(IntNeuros)
	updateNetwork(Theta)
	updateNetwork(Heads)
	
	# Plot data
	plotInterval = 1.0
	if(fmod(t, plotInterval) < dt):
		print('t=%f'%t)
		
		plt.clf()
#		sp = get_neuron_entry(Stellates, 'sp')
		sp = get_neuron_entry(Stellates, 'Vp')
		plt.subplot(2,2,1)
		plt.plot(sp, 'gs')
#		plt.plot(x, dp, 'o')
#		plt.ylim((0,1))
		plt.ylim((-100,100))
		plt.title('Stellate memb. pot.')
		
		plt.subplot(2,2,2)
		sp = get_neuron_entry(Stellates, 'sp')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Stellate syn. pot.')
		
		plt.subplot(2,2,3)
		sp = get_neuron_entry(IntNeuros, 'Vp')
		plt.plot(sp, 'gs')
		plt.ylim((-100,100))
		plt.title('Interneuron memb. pot.')
		
		plt.subplot(2,2,4)
		sp = get_neuron_entry(Heads, 'sp')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Head cells')
#		sp = get_neuron_entry(IntNeuros, 'sp')
#		plt.plot(sp, 'bo')
#		plt.ylim((0,1))
#		plt.title('Interneuron syn. pot.')
		
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
