import sys
sys.path.append('../')

from math import *
import numpy as np
from nlab import *
import random
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
dt = 0.01
nStell = 32		# Number of stellate cells
nInter = 4		# Number of interneurons
nTheta = 0
nHead = 1		# Number of head cells
nHD2S = 1		# Number of head to stellate connections per stellate cell
nSMperS = 2		# Maximal number of submodules per stellate
#nSPerIN = 15	# Number of stellates per submodule

# Connection strengths
sI = 1.5		# Direct current to stellate cells
hdI = 0#1.5		# Direct current to HD cells
thI = 0			# DC to theta oscillator
inI = 1.0
th2s = 0		# theta oscillation input to stellates
s2in = 5.0		# Stellate to interneuron
in2s = 20.0		# Interneuron to stellate
hd2s = 0#10.0		# HD cell to stellate
tMax = 2000





print('Initializing...')



# Create stellate cells
stell = np.array([Neuron_Stel() for _ in range(0,nStell)])
for i in range(0,nStell):
#	stell[i].VP += random.random()*20
#	stell[i].I = sI * (1 + float(i*i)/(nStell*nStell))
	stell[i].I = sI * (1 + i*1.0/nStell)
#	stell[i].I = sI

# Create interneurons and link to one another
inter = np.array([Neuron_Inter() for _ in range(0,nInter)])
connect_many_to_many(inter, inter, in2s)
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

# Distance-determined connections from stellates to interneurons
connRad = (3*nStell)/(4*nInter)
for i in range(nInter):
	k = int(nStell*(i/float(nInter)))
	ind = range(k-connRad, k+connRad)
	for j in ind:
		s2smInd[j].append(i)
	connect_many_to_one(stell[ind], inter[i], s2in)
	connect_one_to_many(inter[i], stell[ind], in2s)

# Create head cells
#head = np.array([Neuron() for _ in range(0,nHead)])
head = np.array([Neuron_Stel() for _ in range(0,nHead)])
for i in range(nHead):
	head[i].I = hdI #*(1 + i/float(nHead))
#	head[i].VP += random.random()*30

# Link HD cells to stellates
s2hdInd = [[] for _ in range(nStell)]	# For each stellate, list of all HDs to which it is connected
for i in range(0,nStell):
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


# Print submodule indices
print('Stellate -> submodule mapping:')
for i in range(nStell):
	print(s2smInd[i])


print('Initialization finished.')
print('Running simulation...')



t = 0.0
m = 0
plotInd = 0



plotLive = False
plotEEG = False
plotFiring = False
plotStelS = False
plotFreq = True

fireHist = [[] for _ in range(nStell+nInter)]
tFireHist = [[] for _ in range(nStell+nInter)]
if plotLive:
	plt.figure()
	plt.ion()
elif plotEEG:
	tHist = []
	sVHist = [[] for _ in range(nStell)]
	sSHist = [[] for _ in range(nStell)]
	hdVHist = [[] for _ in range(nHead)]
	intVHist = [[] for _ in range(nInter)]
	thVHist = [[] for _ in range(nTheta)]

# MAIN TIMELOOP
while(t < tMax):
	t = t+dt
	m = m+1
	
#	# Update HD cells
#	direction = fmod(2*pi*(1 + t/1000), 2*pi)		# Direction in radians, from x-axis
#	for i in range(0, nHead):
#		b = direction/(2*pi)
#		c = 0.25
#		head[i].I = hdI*max(0.0, 1 - 1/c*fabs(fmod(1.0 + float(i)/float(nHead)+c-b, 1)-c))
	
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
	stepNetwork(theta, t, dt)
	stepNetwork(stell, t, dt)
	stepNetwork(inter, t, 2*dt)
	stepNetwork(head, t, dt)
	
	# Check if stellates have fired
	for i in range(nStell):
#		if stell[i].V>0 and stell[i].VP < 50 and stell[i].V > 50:
		if stell[i].isFiring:
			fireHist[i].append(stell[i].V)	
			tFireHist[i].append(t)
	for i in range(nInter):
		if inter[i].isFiring:
			tFireHist[i+nStell].append(t)
	
	# Move to next timestep
	updateNetwork(theta)
	updateNetwork(stell)
	updateNetwork(inter)
	updateNetwork(head)
	
	
	# Plot data live?
	plotInterval = 10.0
	if plotLive and fmod(t, plotInterval) < dt:
#		print('t=%f'%t)
#		print('theta:sp = %f'%theta[0].sp)
	
		plt.clf()
		sp = get_neuron_item(stell, 'V')
		plt.subplot2grid((2,3), (0,0))
		plt.plot(sp, 'gs');		plt.ylim((-70,-10));		plt.title('Stellate memb. pot.')
	
		sp = get_neuron_item(stell, 's')
		plt.subplot2grid((2,3), (0,1))
		plt.plot(sp, 'bo');		plt.ylim((0,1));		plt.title('Stellate syn. pot.')
	
		plt.subplot2grid((2,3), (0,2), rowspan=2)
		plt.plot(tHist, hist, 'rs');		plt.xlim((t-400,t));		plt.title('Firing history')
	
		plt.subplot2grid((2,3), (1,0))
		sp = get_neuron_item(inter, 'V')
		plt.plot(sp, 'gs');		plt.ylim((-70,-10));		plt.title('Interneuron memb. pot.')
	
		plt.subplot2grid((2,3), (1,1))
		sp = get_neuron_item(head, 'V')
		plt.plot(sp, 'gs');		plt.ylim((-70,50));		plt.title('HD memb. pot.')
	
		plt.draw()
		#plt.savefig('./fig/%d.png'%plotInd)
		++plotInd
	
	# Store data for plotting
	saveDataInterval = 0.2
	if plotEEG and fmod(t, saveDataInterval)<dt:
		tHist.append(t)
		for i in range(nStell):
			sVHist[i].append(stell[i].V)
			sSHist[i].append(stell[i].s)
		for i in range(nInter):
			intVHist[i].append(inter[i].V)
		for i in range(nHead):
			hdVHist[i].append(head[i].V)
		for i in range(nTheta):
			thVHist[i].append(theta[i].V)


print('Simulation finished')
print('Plotting data...')



# Compute frequencies
freq = [[] for _ in range(len(tFireHist))]
for i in range(len(tFireHist)):
	nFire = max(len(tFireHist[i])-1, 0)
	freq[i] = np.zeros(nFire)
	for j in range(nFire):
		freq[i][j] = 1000/(tFireHist[i][j+1] - tFireHist[i][j])


if plotEEG or plotFreq:
	nPlot = 0
	if plotEEG:
		nPlot += nStell+nInter+nHead+nTheta
	if plotStelS:
		nPlot += nStell
	if plotFreq:
		nPlot += nStell+nInter
	
	nPlotX = 2
	nPlotY = ceil(float(nPlot)/nPlotX)
	
	plt.figure()
	plt.ion()
	ctr = 1
	for i in range(nStell):
		if plotEEG:
			plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
			plt.plot(tHist,sVHist[i], 'r');		plt.xlim((0, tMax))
	#		for j in range(len(fireHist[i])):
	#			f = fireHist[i][j]
	#			plt.plot(f[0], f[1], 'rs')
	#			# Print instantaneous frequency
	#			if j > 1:
	#				plt.annotate("%d"%(1000.0/(f[0]-fireHist[i][j-1][0])), f, ha='left')
			plt.ylabel("SM:%s. HD:%s"%(stringify(s2smInd[i]), stringify(s2hdInd[i])), rotation="horizontal")
			if plotStelS:
				plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
				plt.plot(tHist,sSHist[i], 'g');			plt.xlim((0, tMax))
		if plotFreq:
			plt.subplot(nPlotY,nPlotX,ctr);		ctr += 1
			if len(freq[i])>0: 
				plt.plot(tFireHist[i][1:], freq[i], 'ro-')
				plt.xlim((0, t))
				plt.ylim((min(freq[i])-1, max(freq[i])+1))
			plt.ylabel('hz, S%d'%i)
	for i in range(nInter):
		if plotEEG:
			plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
			plt.plot(tHist,intVHist[i],'b');	plt.xlim((0, tMax))
			plt.ylabel('hz, IN%d'%i)
		if plotFreq:
			plt.subplot(nPlotY,nPlotX,ctr);		ctr += 1
			if len(freq[i+nStell])>0: 
				plt.plot(tFireHist[i+nStell][1:], freq[i+nStell], 'bo-')
				plt.xlim((0, t))
				plt.ylim((min(freq[i+nStell])-1, max(freq[i+nStell])+1))
			plt.ylabel('Freq')
	if plotEEG:
		for i in range(nHead):
			plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
			plt.plot(tHist,hdVHist[i],'g');		plt.xlim((0, tMax))
		for i in range(nTheta):
			plt.subplot(nPlotY,nPlotX,ctr);	ctr += 1
			plt.plot(tHist, thVHist[i], 'c');		plt.xlim((0, tMax))
	
	plt.draw()


if plotFiring:
	ctr = 0
	plt.figure()
	plt.ion()
	for i in range(nStell):
		plot(tFireHist[i], ctr*np.ones(tFireHist[i].size), 'r.')
		ctr += 1
	plt.ylim((-1, ctr))
