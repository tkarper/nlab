import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt

NStrip = 2**7
dt = 0.005

alpha = 1.0
l = 5
W0 = -1.5
I = 3.0

Head  = array([Neuron_IF() for _ in range(0,2)])			# The Head cells
Strip = array([Neuron_HH() for _ in range(0,2*NStrip)])		# The Strip cells
Right = Strip[0:NStrip]	
Left  = Strip[NStrip:2*NStrip]

for i in range(0,2*NStrip):
	Strip[i].I = I

connect_one_to_many(Head[0], Right, alpha)	# Connecting the Right neurons to the first Head cell
connect_one_to_many(Head[1], Left, alpha)   # Connecting the Left neurons to the second Head cell

M = strip_matrix_OneBump(NStrip, l, W0) 	# Create connection matrix for the Strip cell network
connect_with_matrix(Strip, Strip, M) 	# Apply the connection matrix to the Stripcell network


# INITIAL CONDITIONS
for i in range(NStrip/2, NStrip/2+l):
	Right[i].Vp = 50.0
	Left[i].Vp = 50.0
	Right[i].sp = 0.9
	Left[i].sp  = 0.9
Head[0].sp = 1.0
Head[1].sp = 0.0


t= 0
m = 0
d = 0
plt.ion()

# MAIN TIMELOOP
while(t < 100):
	t= t+dt
	stepNetwork(Strip, dt)  # Perform time-step in Strip cells
	updateNetwork(Strip)	# Update Neuron.sp = Neuron.s
	ll = get_spike_rates(Left)
	rr = get_spike_rates(Right)
	
	if(m%20 == 0):
		plt.clf()
		plt.plot(ll)
		plt.plot(rr)
		plt.title('t=%f'%t)
#		plt.ylim((-80, 50))
		plt.draw()
		#plt.savefig('./fig/%d.png'%d)
		d= d+1
		
#	if(m==80):
#		Head[0].sp = 1.0
#	if(m==1000):
#		Head[1].sp = 1.0
#		Head[0].sp = 0.0
	
	m= m+1
