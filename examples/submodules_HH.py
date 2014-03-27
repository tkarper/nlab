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
I = 1.0			# External exitatory input
s2in = 1.0		# Stellate to interneuron
in2s = -1.0		# Interneuron to stellate
hd2s = 1.0		# HD cell to stellate

# TIF constants: depression decay rate and input rate
cvar.Neuron_TIF_C1 = 0.05
cvar.Neuron_TIF_C2 = 0.1



print('Initializing...')


# Create stellate cells
Stellates = np.array([Neuron_HH() for _ in range(0,nStell)])
for i in range(0,nStell):
	Stellates[i].I = I
cvar.Neuron_HH_alpha = 2
cvar.Neuron_HH_beta = 2
	
# Create interneurons
nIntNeuro = 10	# Number of interneurons
IntNeuros = np.array([Neuron_IF() for _ in range(0,nIntNeuro)])

# Link interneurons to stellate cells and to one another
connect_many_to_many(IntNeuros, IntNeuros, in2s)
nSPerIN = 30	# Number of stellates per submodule
for i in range(0,nIntNeuro):
	submodInd = random.sample(xrange(nStell), nSPerIN)
	submod = Stellates[submodInd]
	connect_many_to_one(submod, IntNeuros[i], s2in)
	connect_one_to_many(IntNeuros[i], submod, in2s)

# Create head cells
nHead = 3		# Number of head cells
Heads = np.array([Neuron_IF() for _ in range(0,nHead)])

# Link HD cells to stellates
Nh2sConns = 3	# Number of head to stellate connections per stellate cell
for i in range(0,nStell):
#	headInd = random.sample(xrange(nHead), Nh2sConns)
#	connect_many_to_one(Heads[headInd], Stellates[i], hd2s)
	for j in range(0,Nh2sConns):
		headInd = random.randrange(nHead)
		Stellates[i].connect(Heads[headInd], hd2s)


print('Initialization finished.')


t = 0.0
m = 0
plotInd = 0
plt.figure()
plt.ion()

direction = 0.0
dirInd = 0.0

# MAIN TIMELOOP
while(t < 100):
	t = t+dt
	m = m+1
	
	# Update HD cells
	direction = 0.0
#	direction = fmod(2*pi + pi/8.0*sin(t/10.0), 2*pi)		# Direction in radians, from x-axis
	for i in range(0, nHead):
		b = direction/(2*pi)
		c = 0.1
		Heads[i].sp = max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
#		Heads[dirInd].sp = 0.0		# Turn off previous HD cell
#		dirInd = round(math.fmod(direction/(2*math.pi), 1.0)*nHead - 0.5)
#		Heads[dirInd].sp = 1.0		# Turn on HD cell
	
	# Update neural network
	stepNetwork(Stellates, t, dt)
	stepNetwork(IntNeuros, t, dt)
	updateNetwork(Stellates)
	updateNetwork(IntNeuros)
	
	# Plot data
	if(m%20 == 0):
		plt.clf()
#		sp = get_neuron_entry(Stellates, 'sp')
		sp = get_neuron_entry(Stellates, 'Vp')
		plt.subplot(2,2,1)
		plt.plot(sp, 'gs')
#		plt.plot(x, dp, 'o')
#		plt.ylim((0,1))
		plt.ylim((-100,100))
		plt.title('t=%f'%t)
		
		plt.subplot(2,2,2)
		sp = get_neuron_entry(Stellates, 'sp')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Stellate s')
		
		plt.subplot(2,2,3)
		sp = get_neuron_entry(Heads, 'sp')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Head cells')
		
		plt.subplot(2,2,4)
		sp = get_neuron_entry(IntNeuros, 'sp')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Interneurons')
		
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
