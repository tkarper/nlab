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
sI = 0.0		# Dirrect current to stellate cells
hdI = 1.5		# Direct current to HD cells
th2s = 0.0		# theta oscillation input to stellates
s2in = 2.0		# Stellate to interneuron
in2s = 5.0		# Interneuron to stellate
hd2s = 4.0		# HD cell to stellate


def get_neuron_item(neurArr, varname):
	arr = np.zeros(neurArr.shape)
	Type = type(neurArr[0])
	for i  in range(neurArr.size):
		arr[i] = Type.__swig_getmethods__[varname](neurArr[i])
	return arr

def set_neuron_item(neurArr, varname, newVal):
	Type = type(neurArr[0])
	for i  in range(neurArr.size):
		Type.__swig_setmethods__[varname](neurArr[i], newVal)



print('Initializing...')



# Create stellate cells
stell = np.array([Neuron_Stel() for _ in range(0,nStell)])
for i in range(0,nStell):
	stell[i].VP += random.random()*10
	stell[i].I = sI

# Create interneurons and link to one another
nInter = int(0.05*nStell)		# Number of interneurons
inter = np.array([Neuron_IntN() for _ in range(0,nInter)])
connect_many_to_many(inter, inter, in2s)

# Link interneurons to stellate cells
nSPerIN = int(0.2*nStell)	# Number of stellates per submodule
# Uniformly random sample
for i in range(0,nInter):
	submodInd = random.sample(xrange(nStell), nSPerIN)
	submod = stell[submodInd]
	connect_many_to_one(submod, inter[i], s2in)
	connect_one_to_many(inter[i], submod, in2s)
### Lateral normal distribution (with replacement!!!)
#distr = np.random.normal(0.0, 0.2, (nStell, nInter))
#distr = distr
#for i in range(0,nInter):
#	# Convert normally distributed numbers to indices in [0, nStell)
#	submodInd = np.round(i*nStell/float(nInter) + nStell*distr[...,i]).astype(int)
#	# Compute indices modulo nStell, avoiding negative index values
#	submodInd = np.mod(nStell + np.mod(submodInd, nStell), nStell)
#	print(np.sort(submodInd))
#	# Extract the submodule and create connections to and from interneuron i
#	submod = stell[submodInd]
#	connect_many_to_one(submod, inter[i], s2in)
#	connect_one_to_many(inter[i], submod, in2s)

# Create head cells
nHead = 2		# Number of head cells
#head = np.array([Neuron() for _ in range(0,nHead)])
head = np.array([Neuron_Stel() for _ in range(0,nHead)])
head[0].I = hdI

# Link HD cells to stellates
nHD2S = 1	# Number of head to stellate connections per stellate cell
for i in range(0,nStell):
	headInd = random.sample(xrange(nHead), nHD2S)
	connect_many_to_one(head[headInd], stell[i], hd2s)
#	for j in range(0,Nh2sConns):
#		headInd = random.randrange(nHead)
#		stell[i].connect(head[headInd], hd2s)

# Create theta oscillator and link to stellates
theta = np.array([Neuron_Osc(50, 3, 1)])	# Period, duration, strength
connect_one_to_many(theta[0], stell, th2s)
#for i in range(0,nStell):
#	stell[i].I = I


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
#		head[i].s = max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
	
	# Update neural network
	stepNetwork(theta, t, dt)
	stepNetwork(stell, t, dt)
	stepNetwork(inter, t, dt)
	stepNetwork(head, t, dt)
	
	# Check if stellates have fired
	for i in range(0, nStell):
		if stell[i].V>0 and stell[i].V < stell[i].VP:
				tHist.append(t)
				hist.append(i)
	
	# Move to next timestep
	updateNetwork(stell)
	updateNetwork(inter)
	updateNetwork(theta)
	updateNetwork(head)
	
	# Turn off external input for an initial phase
	if (t < 10):
		theta[0].sp = 0
		for i in range(0,nHead):
			head[i].sp = 0
	
	# Plot data
	plotInterval = 10.0
	if(fmod(t, plotInterval) < dt):
#		print('t=%f'%t)
#		print('theta:sp = %f'%theta[0].sp)
		
		plt.clf()
		sp = get_neuron_item(stell, 'V')
#		plt.subplot(2,3,1)
		plt.subplot2grid((2,3), (0,0))
		plt.plot(sp, 'gs')
		plt.ylim((-70,-10))
		plt.title('Stellate memb. pot.')
		
#		plt.subplot(2,3,2)
		sp = get_neuron_item(stell, 's')
		plt.subplot2grid((2,3), (0,1))
		plt.plot(sp, 'bo')
		plt.ylim((0,1))
		plt.title('Stellate syn. pot.')
		
#		plt.subplot(2,3,5)
#		sp = get_neuron_item(head, 's')
#		plt.plot(sp, 'bo')
#		plt.ylim((0,1))
#		plt.title('Head cells')
		plt.subplot2grid((2,3), (0,2), rowspan=2)
		plt.plot(tHist, hist, 'rs')
		plt.xlim((t-400,t))
		plt.title('Firing history')
		
#		plt.subplot(2,3,4)
		plt.subplot2grid((2,3), (1,0))
		sp = get_neuron_item(inter, 'V')
		plt.plot(sp, 'gs')
		plt.ylim((-70,-10))
		plt.title('Interneuron memb. pot.')
		
		
		plt.subplot2grid((2,3), (1,1))
		sp = get_neuron_item(head, 'V')
		plt.plot(sp, 'gs')
		plt.ylim((-70,50))
		plt.title('HD memb. pot.')
		
		
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
