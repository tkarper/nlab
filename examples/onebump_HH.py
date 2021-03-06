import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt

NStrip = 2**6
dt = 0.005

alpha = 1.0
l = 10
W0 = -1.0
I = 1.0

#Head  = array([Neuron_IF() for _ in range(0,2)])			# The Head cells
Theta = Neuron_Osc(30, 1, 1)						# Theta oscillation
Head = array([Neuron_Osc(5, 2, 1) for _ in range(0,2)])			# The Head cells
InputNeurons = array([Head[0], Head[1], Theta])
Strip = array([Neuron_HH() for _ in range(0,2*NStrip)])		# The Strip cells
Right = Strip[0:NStrip]	
Left  = Strip[NStrip:2*NStrip]

for i in range(0,2*NStrip):
	Strip[i].I = I
cvar.Neuron_HH_VT = 0
cvar.Neuron_HH_alpha = 1
cvar.Neuron_HH_beta = 0.4

connect_one_to_many(Theta, Strip, I)
connect_one_to_many(Head[0], Right, alpha)	# Connecting the Right neurons to the first Head cell
connect_one_to_many(Head[1], Left, alpha)   # Connecting the Left neurons to the second Head cell

M = strip_matrix_OneBump2(NStrip, l, W0) 	# Create connection matrix for the Strip cell network
connect_with_matrix2(Strip, Strip, M) 	# Apply the connection matrix to the Stripcell network


# INITIAL CONDITIONS
for i in range(NStrip/2, NStrip/2+1):
	Right[i].Vp = -0.0
	Left[i].Vp = -0.0
#	Right[i].sp = 1.0
#	Left[i].sp  = 1.0
#Head[0].sp = 1.0
#Head[1].sp = 0.0
Head[0].strength = 1
Head[1].strength = 0


t= 0
m = 0
d = 0
plt.figure()
plt.ion()
plotVarName = "Vp"; yaxis = (-100, 70)
#plotVarName = "np"; yaxis = (-0.1, 1.1)
#plotVarName = "sp"; yaxis = (-0.1, 1.1)
# MAIN TIMELOOP
while(t < 200):
	t= t+dt
	stepNetwork(InputNeurons, t, dt)
	updateNetwork(InputNeurons)
	stepNetwork(Strip, t, dt)  # Perform time-step in Strip cells
	updateNetwork(Strip)	# Update Neuron.sp = Neuron.s
	ll = get_neuron_entry(Left, plotVarName)
	rr = get_neuron_entry(Right, plotVarName)
	
	if(m%40 == 0):
		plt.clf()
		plt.plot(ll, '-o')
		plt.plot(rr, '-s')
		plt.title('t=%f'%t)
		plt.ylim(yaxis)
		plt.draw()
		#plt.savefig('./fig/%d.png'%d)
		d= d+1
		
	# if(m==80):
	# 	Head[0].sp = 1.0
	if(abs(t-50) < dt):
		Head[0].strength = 0
		Head[1].strength = 1
	
	m= m+1
