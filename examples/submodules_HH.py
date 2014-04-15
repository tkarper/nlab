import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
import time
import matplotlib.pyplot as plt


# UTILITY FUNCTIONS
# Get all elements of an array of SWIG objects
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





# Discretization parameters
dt = 0.02
nStell = 32		# Number of stellate cells
nInter = 1		# Number of interneurons
nTheta = 0
nHead = 0		# Number of head cells
nHD2S = 3		# Number of head to stellate connections per stellate cell
nSMperS = 2		# Maximal number of submodules per stellate
#nSPerIN = 15	# Number of stellates per submodule

# Connection strengths
sI = 2.0		# Direct current to stellate cells
inI = 0.4
hdI = 0			# Direct current to HD cells
thI = 0			# DC to theta oscillator
th2s = 0		# theta oscillation input to stellates
s2in = 2.0		# Stellate to interneuron
in2s = 2.0		# Interneuron to stellate
in2in= in2s
hd2s = 0.2		# HD cell to stellate
tMax = 10000

# Seed random number generator
randSeed = 273#int(time.time()) % 10000
#randSeed = 4145
print('Random seed: %d'%randSeed)
random.seed(randSeed)





print('Initializing...')



# Create stellate cells
stell = np.array([Neuron_Traub() for _ in range(0,nStell)])
#stell = np.array([Neuron_Stel() for _ in range(0,nStell)])
for i in range(0,nStell):
#	stell[i].VP += random.random()*20
#	stell[i].I = sI * (1 + float(i*i)/(nStell*nStell))
#	stell[i].I = sI * (1 + i*0.1/nStell)
#	stell[i].I = sI * (1 + 0.2*random.random())
	stell[i].I = sI * random.random()
#	stell[i].I = sI
	stell[i].gM = 1

# Create interneurons and link to one another
inter = np.array([Neuron_Traub_IN() for _ in range(0,nInter)])
#inter = np.array([Neuron_Inter() for _ in range(0,nInter)])
connect_many_to_many(inter, inter, in2in)
for i in range(nInter):
	inter[i].I = inI

# Create submodules
s2smInd = [[] for _ in range(nStell)]	# For each stellate, list of all submodules to which it belongs
#submodInd = [[] for _ in range(nInter)]
## Uniformly random sample
#for i in range(nInter):
#	submodInd[i] = random.sample(xrange(nStell), nSPerIN)
# Lateral normal distribution (with replacement!!!)
#distr = np.random.normal(0.0, 1.0/nInter, (nStell, nInter))
#for i in range(nInter):
#	# Convert normally distributed numbers to indices in [0, nStell)
#	submodInd[i] = np.round(i*nStell/float(nInter) + nStell*distr[...,i]).astype(int)
#	# Compute indices modulo nStell, avoiding negative index values
#	submodInd[i] = np.unique(np.mod(nStell + np.mod(submodInd[i], nStell), nStell))
#
## Link interneurons to stellate cells
#for i in range(nInter):
#	# Extract the submodule and create connections to and from interneuron i
#	submod = stell[submodInd[i]]
#	connect_many_to_one(submod, inter[i], s2in)
#	connect_one_to_many(inter[i], submod, in2s)
#	for ind in submodInd[i]:
#		s2smInd[ind].append(i)

## Normally distributed random connections from stellates to interneurons
#distr = np.random.normal(0.0, 0.5/nInter, (nStell, nSMperS))
#for i in range(nStell):
#	# Convert normally distributed numbers to indices in [0, nInter)
#	ind = np.round(nInter*(i/float(nStell) + distr[i,...])).astype(int)
#	# Compute indices modulo nInter and remove duplicates
#	ind = np.unique(np.mod(np.mod(nInter+ind,nInter), nInter))
#	# Link interneurons to stellate cells
#	connect_many_to_one(inter[ind], stell[i], in2s)
#	connect_one_to_many(stell[i], inter[ind], s2in)
#	s2smInd[i] = ind

## Per-stellate based random connections
#for i in range(nStell):
#	ind = np.random.randint(0, nInter, nSMperS)
#	for j in ind:
#		s2smInd[i].append(j)
#	connect_one_to_many(stell[i], inter[ind], s2in)
#	connect_many_to_one(inter[ind], stell[i], in2s)

# Distance-determined connections from stellates to interneurons
connRad = nStell/(2*nInter)
#connRad = nStell/nInter
#connRad = (3*nStell)/(2*nInter)
for i in range(nInter):
	k = int(nStell*(i/float(nInter)))
	ind = range(k-connRad, k+connRad)
	for j in range(len(ind)):
		ind[j] = ind[j] % nStell
		s2smInd[ind[j]].append(i)
	connect_many_to_one(stell[ind], inter[i], s2in)
	connect_one_to_many(inter[i], stell[ind], in2s)

# Create head cells
head = np.array([Neuron_Traub() for _ in range(nHead)])
for i in range(nHead):
	head[i].I = hdI #*(1 + i/float(nHead))
#	head[i].VP += random.random()*30

# Link HD cells to stellates
s2hdInd = [[] for _ in range(nStell)]	# For each stellate, list of all HDs to which it is connected
if nHead > 0:
	for i in range(nStell):
		headInd = random.sample(xrange(nHead), nHD2S)
		connect_many_to_one(head[headInd], stell[i], hd2s)
		s2hdInd[i] = headInd
	#	stell[i].connect(stell[random.randint(0,nStell-1)], s2in)
	#	for j in range(0,Nh2sConns):
	#		headInd = random.randrange(nHead)
	#		stell[i].connect(head[headInd], hd2s)

# Create theta oscillator and link to stellates
# Create theta oscillator
theta = np.array([Neuron_Stel() for _ in range(nTheta)])
for th in theta:
	th.I = thI
connect_many_to_many(theta, stell, th2s)
connect_many_to_many(theta, head, th2s)





allNeurons = np.concatenate((stell, inter, head, theta))
nNeuro = allNeurons.size

# Print submodule indices
print('Stellate -> submodule mapping:')
for i in range(nStell):
	print(s2smInd[i])






print('Initialization finished.')
print('Running simulation...')




plotLive = 	0
plotEEG = 	0
plotFiring=	0
plotS = 	0
plotFreq = 	1

tFireHist = [[] for _ in range(nNeuro)]
if plotLive:
	plt.figure()
	plt.ion()
elif plotEEG:
	tHist = []
	VHist = [[] for _ in range(nNeuro)]
	sHist = [[] for _ in range(nNeuro)]

# MAIN TIMELOOP
t = 0.0
m = 0
plotInd = 0
while(t < tMax):
	t = t+dt
	m = m+1
	
	# Update HD cells
	direction = fmod(2*pi*(1 + t/30000), 2*pi)		# Direction in radians, from x-axis
	for i in range(0, nHead):
		b = direction/(2*pi)
		c = 0.2
		head[i].I = 0.4 + hdI*max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
	
#	# Turn off external input for an initial phase
#	if (t < 10):
##		theta[0].sp = 0
#		for i in range(0,nHead):
#			head[i].sp = 0

#	if t>1000 and t<2000:
#	stell[0].I = sI * (1+int(t/1000))
#	stell[0].I = sI * (1+int(t/1000) + float(0*0)/(nStell*nStell))
#	elif t>=2000:
#		stell[0].I = sI * (1 + float(0*0)/(nStell*nStell))


	# Update neural network
	stepNetwork(allNeurons, t, dt)	
	# Move to next timestep
	updateNetwork(allNeurons)
	
	
	# Check if stellates have fired
	for i in range(nNeuro):
		if allNeurons[i].isFiring:
			tFireHist[i].append(t)
	
	# Plot data live?
	plotInterval = 10.0
	if plotLive and fmod(t, plotInterval) < dt:
		plt.clf()
		sp = get_neuron_item(stell, 'V')
		plt.subplot2grid((2,3), (0,0));		plt.plot(sp, 'gs');		plt.ylim((-70,-10));		plt.title('Stellate memb. pot.')
		sp = get_neuron_item(stell, 's')
		plt.subplot2grid((2,3), (0,1));		plt.plot(sp, 'bo');		plt.ylim((0,1));		plt.title('Stellate syn. pot.')
		plt.subplot2grid((2,3), (0,2), rowspan=2)
		plt.plot(tHist, hist, 'rs');		plt.xlim((t-400,t));		plt.title('Firing history')
		plt.subplot2grid((2,3), (1,0))
		sp = get_neuron_item(inter, 'V');		plt.plot(sp, 'gs');		plt.ylim((-70,-10));		plt.title('Interneuron memb. pot.')
		plt.subplot2grid((2,3), (1,1))
		sp = get_neuron_item(head, 'V');		plt.plot(sp, 'gs');		plt.ylim((-70,50));		plt.title('HD memb. pot.')
		plt.draw()
		++plotInd
	
	# Store data for plotting
	saveDataInterval = 0.1
	if plotEEG and fmod(t, saveDataInterval)<dt:
		tHist.append(t)
		for i in range(nNeuro):
			VHist[i].append(allNeurons[i].V)
			sHist[i].append(allNeurons[i].s)




print('Simulation finished')
print('Plotting data...')



# Compute frequencies
smoothFreq = 1		# When computing frequencies, take average over subsequent 'smoothFreq' number of peaks
freq = [[] for _ in range(len(tFireHist))]
for i in range(len(tFireHist)):
	nFire = max(len(tFireHist[i])-1, 0)
	freq[i] = np.zeros(nFire)
	for j in range(nFire):
		k = min(j+smoothFreq, nFire) - j
		freq[i][j] = k*1000/(tFireHist[i][j+k] - tFireHist[i][j])



plotColors = ['r' for _ in range(nStell)] + ['b' for _ in range(nInter)] + ['g' for _ in range(nHead)]
if plotEEG or plotFreq:
	nPlot = 0
	if plotEEG:
		nPlot += nStell+nInter+nHead+nTheta
	if plotS:
		nPlot += nStell+nInter+nHead
	if plotFreq:
		nPlot += nStell+nInter+nHead
	
	nPlotX = 3
	nPlotY = ceil(float(nPlot)/nPlotX)
	
	plt.figure()
	plt.ion()
	ctr = 1
	for i in range(nNeuro):
		if plotEEG:
			plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
			plt.plot(tHist,VHist[i], plotColors[i]);		plt.xlim((0, tMax))
#			plt.ylabel("SM:%s. HD:%s"%(stringify(s2smInd[i]), stringify(s2hdInd[i])), rotation="horizontal")
			if plotS:
				plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
				plt.plot(tHist, sHist[i], 'g');			plt.xlim((0, tMax))
		if plotFreq:
			plt.subplot(nPlotY,nPlotX,ctr);		ctr += 1
			if len(freq[i])>0: 
				plt.plot(tFireHist[i][1:], freq[i], plotColors[i]+'o')
				plt.xlim((0, t))
				plt.ylim((min(freq[i])-1, max(freq[i])+1))
			if i<nStell:
				label = "SM:%s"%stringify(s2smInd[i])
				if nHead > 0:
					label += ". HD:%s"%stringify(s2hdInd[i])
				plt.ylabel(label, rotation="horizontal")
	plt.draw()


if plotFiring:
	ctr = 0
	plt.figure()
	plt.ion()
	for i in range(nStell):
		plot(tFireHist[i], ctr*np.ones(tFireHist[i].size), 'r.')
		ctr += 1
	plt.ylim((-1, ctr))
