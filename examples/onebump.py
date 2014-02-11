import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt

NStrip = 16
dt = 0.1

alpha = 0.5
l  = 3
W0 = -1.0

Head  = array([Neuron_IF() for _ in range(0,2)])			# The Head cells
Strip = array([Neuron_IF() for _ in range(0,2*NStrip)])		# The Strip cells
Right = Strip[0:NStrip]	
Left  = Strip[NStrip:2*NStrip]


connect_one_to_many(Head[0], Right, alpha)	# Connecting the Right neurons to the first Head cell
connect_one_to_many(Head[1], Left, alpha)   # Connecting the Left neurons to the second Head cell

M = stripCell_onebump(Right, Left, l,W0) 	# Create connection matrix for the Strip cell network
connect_with_matrix(Strip, Strip, M) 	# Apply the connection matrix to the Stripcell network


#INITIAL CONDITIONS
Right[NStrip/2].sp = 1.0
Left[NStrip/2].sp  = 1.0


t= 0
m = 0
d = 0
plt.ion()

# while(t < 20):
# 	t= t+dt
# 	m= m+1
# 	stepNetwork(Head, dt)	# Perform time-step in Head cells
# 	stepNetwork(Strip, dt)  # Perform time-step in Strip cells
# 	updateNetwork(Head)		# Update Neuron.sp = Neuron.s
# 	updateNetwork(Strip)	# Update Neuron.sp = Neuron.s
# 	ll = get_spike_rates(Left)
# 	rr = get_spike_rates(Right)
# 	
# 	if(m%1 == 0):
# 		plt.clf()
# 		plt.plot(ll)
# 		plt.plot(rr)
# 		plt.savefig('./fig/%d.png'%d)
# 		d= d+1