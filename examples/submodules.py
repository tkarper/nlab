import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

# Discretization parameters
nStell = 1000	# Number of stellate cells
dt = 0.01

# Connection strengths
I = 0.5			# External exitatory input
s2in = 1.0		# Stellate to interneuron
in2s = -1.0		# Interneuron to stellate
hd2s = 1.0		# HD cell to stellate

# Other constants
cvar.Neuron_TIF_C1 = 0.05
cvar.Neuron_TIF_C2 = 0.1



print('Initializing...')


# Create stellate cells
Stellates = np.array([Neuron_TIF() for _ in range(0,nStell)])
for i in range(0,nStell):
	Stellates[i].I = I
#	Stellates[i].sp = random.random()*0.1		# Initialize randomly
	
# Create interneurons
nIntNeuro = 10	# Number of interneurons
IntNeuro = np.array([Neuron_IF() for _ in range(0,nIntNeuro)])

# Link interneurons to stellate cells and to one another
connect_many_to_many(IntNeuro, IntNeuro, in2s)
nSPerIN = 700	# Number of stellates per submodule
for i in range(0,nIntNeuro):
	submodInd = random.sample(xrange(nStell), nSPerIN)
	submod = Stellates[submodInd]
	connect_many_to_one(submod, IntNeuro[i], s2in)
	connect_one_to_many(IntNeuro[i], submod, in2s)

# Create head cells
nHead = 32		# Number of head cells
Heads = np.array([Neuron_IF() for _ in range(0,nHead)])

# Link HD cells to stellates
Nh2sConns = 3	# Number of head to stellate connections per stellate cell
for i in range(0,nStell):
	stell = Stellates[i]
	for j in range(0,Nh2sConns):
		headInd = random.randrange(nHead)
		stell.connect(Heads[headInd], hd2s)


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
	stepNetwork(IntNeuro, t, dt)
	updateNetwork(Stellates)
	updateNetwork(IntNeuro)
	
	# Plot data
	if(m%20 == 0):
		sp = get_neuron_entry(IntNeuro, 'sp')
		print(sp)
		
		plt.clf()
		sp = get_neuron_entry(Stellates, 'sp')
		x = np.linspace(0, 1, nStell)
		plt.subplot(1,2,1)
		plt.plot(x, sp, 'gs')
#		plt.plot(x, dp, 'o')
		plt.ylim((0,1))
		plt.title('t=%f'%t)
		
		plt.subplot(1,2,2)
		sp = get_neuron_entry(Heads, 'sp')
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
