import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import matplotlib.pyplot as plt

# Discretization parameters
dt = 0.01
nStell = 24		# Number of stellate cells
nInter = 1		# Number of interneurons
nSPerIN = 15	# Number of stellates per submodule
nHead = 4		# Number of head cells
nHD2S = 1		# Number of head to stellate connections per stellate cell

# Connection strengths
sI = 1.0		# Dirrect current to stellate cells
hdI = 10.0		# Direct current to HD cells
th2s = 0.0		# theta oscillation input to stellates
s2in = 12.0		# Stellate to interneuron
in2s = 15.0		# Interneuron to stellate
hd2s = 15.0		# HD cell to stellate
tMax = 2000


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

# Stringify a list of integers
def stringify(x):
	label = ""
	separator = ""
	for j in x:
		label = "%s%s%i"%(label, separator, j)
		separator = ","
	return label



print('Initializing...')



# Create stellate cells
stell = np.array([Neuron_Stel() for _ in range(0,nStell)])
for i in range(0,nStell):
	stell[i].VP += random.random()*20
	stell[i].I = sI

# Create interneurons and link to one another
inter = np.array([Neuron_IntN() for _ in range(0,nInter)])
connect_many_to_many(inter, inter, in2s)

# Create submodules
s2smInd = [[] for _ in range(nStell)]	# For each stellate, list of all submodules to which it belongs
submodInd = [[] for _ in range(nInter)]
# Uniformly random sample
for i in range(nInter):
	submodInd[i] = random.sample(xrange(nStell), nSPerIN)
# Lateral normal distribution (with replacement!!!)
#distr = np.random.normal(0.0, 1.0/nInter, (nStell, nInter))
#distr = distr
#for i in range(nInter):
#	# Convert normally distributed numbers to indices in [0, nStell)
#	submodInd[i] = np.round(i*nStell/float(nInter) + nStell*distr[...,i]).astype(int)
#	# Compute indices modulo nStell, avoiding negative index values
#	submodInd[i] = np.unique(np.mod(nStell + np.mod(submodInd[i], nStell), nStell))

# Link interneurons to stellate cells
for i in range(nInter):
	# Extract the submodule and create connections to and from interneuron i
	submod = stell[submodInd[i]]
	connect_many_to_one(submod, inter[i], s2in)
	connect_one_to_many(inter[i], submod, in2s)
	for ind in submodInd[i]:
		s2smInd[ind].append(i)

# Create head cells
#head = np.array([Neuron() for _ in range(0,nHead)])
head = np.array([Neuron_Stel() for _ in range(0,nHead)])
#head[0].I = hdI
#head[0].mDepRise = 0

# Link HD cells to stellates
s2hdInd = [[] for _ in range(nStell)]	# For each stellate, list of all HDs to which it is connected
for i in range(0,nStell):
	headInd = random.sample(xrange(nHead), nHD2S)
	connect_many_to_one(head[headInd], stell[i], hd2s)
	s2hdInd[i] = headInd
#	for j in range(0,Nh2sConns):
#		headInd = random.randrange(nHead)
#		stell[i].connect(head[headInd], hd2s)

## Create theta oscillator and link to stellates
#theta = np.array([Neuron_Osc(50, 3, 1)])	# Period, duration, strength
#connect_one_to_many(theta[0], stell, th2s)
##for i in range(0,nStell):
##	stell[i].I = I



print('Initialization finished.')
print('Starting simulation')



t = 0.0
m = 0
plotInd = 0



plotLive = False

fireHist = [[] for _ in range(nStell)]
if plotLive:
	plt.figure()
	plt.ion()
else:
	tHist = []
	sVHist = [[] for _ in range(nStell)]
	sSHist = [[] for _ in range(nStell)]
	hdVHist = [[] for _ in range(nHead)]
	intVHist = [[] for _ in range(nInter)]

# MAIN TIMELOOP
while(t < tMax):
	t = t+dt
	m = m+1
	
	# Update HD cells
#	direction = 0.0
	direction = fmod(2*pi*(1 + t/1000), 2*pi)		# Direction in radians, from x-axis
	for i in range(0, nHead):
		b = direction/(2*pi)
		c = 0.25
		head[i].I = hdI*max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
	
	# Update neural network
#	stepNetwork(theta, t, dt)
	stepNetwork(stell, t, dt)
	stepNetwork(inter, t, dt)
	stepNetwork(head, t, dt)
	
	# Check if stellates have fired
	for i in range(0, nStell):
		if stell[i].V>0 and stell[i].VP < 50 and stell[i].V > 50:
			fireHist[i].append((t,stell[i].V))
	
	# Move to next timestep
	updateNetwork(stell)
	updateNetwork(inter)
#	updateNetwork(theta)
	updateNetwork(head)
	
	# Turn off external input for an initial phase
	if (t < 10):
#		theta[0].sp = 0
		for i in range(0,nHead):
			head[i].sp = 0
	
	# Plot data live?
	plotInterval = 10.0
	if plotLive and fmod(t, plotInterval) < dt:
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
	
	if not plotLive and m%10==0:
		tHist.append(t)
		for i in range(nStell):
			sVHist[i].append(stell[i].V)
			sSHist[i].append(stell[i].s)
		for i in range(nInter):
			intVHist[i].append(inter[i].V)
		for i in range(nHead):
			hdVHist[i].append(head[i].V)


print('Simulation finished')


# PLOT DATA
if not plotLive:
	plotStelS = False
	plt.figure()
	plt.ion()
	if plotStelS:
		nPlot = 2*nStell+nInter+nHead
	else:
		nPlot = nStell+nInter+nHead
	ctr = 1
	for i in range(nStell):
		plt.subplot(nPlot,1,ctr);	ctr += 1
		plt.plot(tHist,sVHist[i], 'r')
		plt.xlim((0, tHist[-1]))
		for j in range(len(fireHist[i])):
			f = fireHist[i][j]
			plt.plot(f[0], f[1], 'rs')
#			# Print instantaneous frequency
#			if j > 1:
#				plt.annotate("%d"%(1000.0/(f[0]-fireHist[i][j-1][0])), f, ha='left')
		plt.ylabel("SM:%s. HD:%s"%(stringify(s2smInd[i]), stringify(s2hdInd[i])), rotation="horizontal")
		if plotStelS:
			plt.subplot(nPlot,1,ctr);	ctr += 1
			plt.plot(tHist,sSHist[i], 'g')
			plt.xlim((0, tHist[-1]))
	for i in range(nInter):
		plt.subplot(nPlot,1,ctr);	ctr += 1
		plt.plot(tHist,intVHist[i],'b')
	for i in range(nHead):
		plt.xlim((0, tHist[-1]))
		plt.subplot(nPlot,1,ctr);	ctr += 1
		plt.plot(tHist,hdVHist[i],'g')
		plt.xlim((0, tHist[-1]))
