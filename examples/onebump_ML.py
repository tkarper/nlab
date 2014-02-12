import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt

NStrip = 2**6
dt = 0.005

alpha = 0.0
l = 5
W0 = -0.0
I = 0.0

Head  = array([Neuron_IF() for _ in range(0,2)])			# The Head cells
Strip = array([Neuron_HH() for _ in range(0,2*NStrip)])		# The Strip cells
Right = Strip[0:NStrip]	
Left  = Strip[NStrip:2*NStrip]

for i in range(0,2*NStrip):
	Strip[i].I = I
cvar.Neuron_HH_VT = 0
cvar.Neuron_HH_beta = 1

connect_one_to_many(Head[0], Right, alpha)	# Connecting the Right neurons to the first Head cell
connect_one_to_many(Head[1], Left, alpha)   # Connecting the Left neurons to the second Head cell

# M = strip_matrix_OneBump(NStrip, l, W0) 	# Create connection matrix for the Strip cell network
# connect_with_matrix(Strip, Strip, M) 	# Apply the connection matrix to the Stripcell network


# INITIAL CONDITIONS
for i in range(NStrip/2, NStrip/2+1):
	Right[i].Vp = -0.0
	Left[i].Vp = -0.0
	Right[i].sp = 0.0
	Left[i].sp  = 0.0
Head[0].sp = 1.0
Head[1].sp = 0.0


t= 0
m = 0
d = 0
plt.ion()
stepNetwork(Strip, dt)  # Perform time-step in Strip cells
updateNetwork(Strip)	# Update Neuron.sp = Neuron.s
ll = get_neuron_entry(Left, "mp")
rr = get_neuron_entry(Left, "np")

	
# MAIN TIMELOOP
while(t < 100):
	t= t+dt
	stepNetwork(Strip, dt)  # Perform time-step in Strip cells
	updateNetwork(Strip)	# Update Neuron.sp = Neuron.s
	ll = get_neuron_entry(Left, "Vp")
	rr = get_neuron_entry(Right, "Vp")
	
	if(m%1 == 0):
		plt.clf()
		plt.plot(ll)
		plt.plot(rr)
		plt.title('t=%f'%t)
		plt.ylim((-80, 50))
		#plt.ylim((0, 1))
		plt.draw()
		#plt.savefig('./fig/%d.png'%d)
		d= d+1
		
	# if(m==80):
	# 	Head[0].sp = 1.0
	# if(m==1000):
	# 	Head[1].sp = 1.0
	# 	Head[0].sp = 0.0
	
	m= m+1
