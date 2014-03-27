import sys
sys.path.append('../')

import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

# Discretization parameters
NStell = 200	# Number of stellate cells
dt = 0.01

# Connection strengths
I = 0.1			# External exitatory input
s2in = 0.5		# Stellate to interneuron
in2s = -2.0		# Interneuron to stellate
hd2s = 1.0		# HD cell to stellate

# Other constants
cvar.Neuron_TIF_C1 = 0.1
cvar.Neuron_TIF_C2 = 0.2



print('Initializing...')


# Create stellate cells
Stellates = np.array([Neuron_TIF() for _ in range(0,NStell)])
for i in range(0,NStell):
	Stellates[i].I = I
	Stellates[i].sp = random.random()*0.1		# Initialize randomly
	
# Create interneurons and link to stellate cells
NIntNeuro = 5	# Number of interneurons
NSPerIN = 100	# Number of stellates per submodule
IntNeuro = np.array([Neuron_TIF() for _ in range(0,NIntNeuro)])
for i in range(0,NIntNeuro):
	submodInd = random.sample(xrange(NStell), NSPerIN)
	submod = Stellates[submodInd]
	print('connecting s2in')
	connect_many_to_many(submod, IntNeuro[i], s2in)
	print('connecting in2s')
	connect_many_to_many(IntNeuro[i], submod, in2s)

# Create head cells
NHead = 16		# Number of head cells
Heads = np.array([Neuron_IF() for _ in range(0,16)])

# Link HD cells to stellates
Nh2sConns = 3	# Number of head to stellate connections per stellate cell
for i in range(0,NStell):
	nStell = Stellates[i]
	for j in range(0,Nh2sConns):
		headInd = random.randrange(NHead)
		nStell.connect(Heads[headInd], hd2s)


print('Initialization finished.')


t = 0
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
	Heads[dirInd].sp = 0.0		# Turn off previous HD cell
	direction = pi/6.0*sin(t/5.0)		# Direction in radians, from x-axis
	dirInd = round(fmod(direction/(2*pi), 1.0)*NHead - 0.5)
	Heads[dirInd].sp = 1.0		# Turn on HD cell
	
	# Update neural network
	stepNetwork(Stellates, t, dt)
	updateNetwork(Stellates)
	
	# Plot data
	if(m%20 == 0):
		plt.clf()
		sp = get_neuron_entry(Stellates, 'sp')
		dp = get_neuron_entry(Stellates, 'dp')
		x = linspace(0, 1, NStell)
		plt.plot(x, sp, 'gs')
#		plt.plot(x, dp, 'o')
		plt.ylim((0,1))
		plt.title('t=%f'%t)
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
