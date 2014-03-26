import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt

NStrip = 2**6
dt = 0.1

alpha = 1
l  = int(ceil(0.1*NStrip))
W0 = -20.0/NStrip
I = 0.0
speed = 1.0
accel = 0.2
velocity = array([speed])


Head  = array([Neuron_IF() for _ in range(0,2)])			# The Head cells
Strip = array([Neuron_IF() for _ in range(0,2*NStrip)])		# The Strip cells
Right = Strip[0:NStrip]	
Left  = Strip[NStrip:2*NStrip]

for i in range(0,2*NStrip):
	Strip[i].I = I

connect_one_to_many(Head[0], Right, alpha)	# Connecting the Right neurons to the first Head cell
connect_one_to_many(Head[1], Left, alpha)   # Connecting the Left neurons to the second Head cell

M = strip_matrix_OneBump2(NStrip, l, W0) 	# Create connection matrix for the Strip cell network
connect_with_matrix2(Strip, Strip, M) 	# Apply the connection matrix to the Stripcell network


# INITIAL CONDITIONS
for i in range(NStrip/2, NStrip/2+l):
	Right[i].sp = 1.0
	Left[i].sp  = 0.0
Head[0].sp = speed


t= 0
m = 0
d = 0
plt.figure()
plt.ion()

# MAIN TIMELOOP
while(t < 100):
	t= t+dt
	m= m+1
	stepNetwork(Strip, t, dt)  # Perform time-step in Strip cells
	updateNetwork(Strip)	# Update Neuron.sp = Neuron.s
	ll = get_spike_rates(Left)
	rr = get_spike_rates(Right)
	
	if(m%1 == 0):
		x = linspace(0, 1, NStrip)
		plt.clf()
		plt.plot(x, ll, '-o')
		plt.plot(x, rr, '-s')
		plt.ylim((0,1))
		plt.title('t=%f'%t)
		plt.draw()
		#plt.savefig('./fig/%d.png'%d)
		d= d+1
	
	# Update velocity randomly
	print("Velocity: %f" % velocity[0])
	dv = (2*random.random() - 1)
	dv = dv/linalg.norm(dv)*accel*dt
	velocity = velocity + dv
	print(velocity[0])
	print("New velocity: %f" % velocity[0])
	Head[0].sp = max(velocity[0], 0.0)
	Head[1].sp = max(-velocity[0], 0)
