import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt
from mayavi import mlab

NStrip = 2**5
dt = 0.1
xGrid, yGrid = mgrid[0:1:1.0/NStrip, 0:1:1.0/NStrip]
xStrip = linspace(0, 1, NStrip)
position = array([0.0,0.0])
posXLim = array([-10.0, 10.0])
posYLim = array([-10.0, 10.0])

alpha = 0.1
l  = 3
W0 = -0.5
I = 1.0
speed = 1.0
accel = 1.0
velocity = array([speed, 0])

# Create Head cells
Head   = array([Neuron_IF() for _ in range(0,4)])			# The Head cells

# Create Strip Cells
Strip  = array([Neuron_IF() for _ in range(0,4*NStrip)])		# The Strip cells
StripX = Strip[2*NStrip:4*NStrip]
StripY = Strip[0:2*NStrip]
Right  = StripX[0:NStrip]	
Left   = StripX[NStrip:2*NStrip]
Up     = StripY[0:NStrip]	
Down   = StripY[NStrip:2*NStrip]

# Create a grid cell module
Grid   = array([Neuron_IF() for _ in range(0,NStrip*NStrip)])

# Connect Head cells to Strip cells
connect_one_to_many(Head[0], Right, alpha)	# Connecting the Right neurons to the first Head cell
connect_one_to_many(Head[1], Left, alpha)   # Connecting the Left neurons to the second Head cell
connect_one_to_many(Head[2], Up, alpha)     # Connecting the Up neurons to the third Head cell
connect_one_to_many(Head[3], Down, alpha)   # Connecting the Down neurons to the fourth Head cell

# Connect the internal structure in strip cells
MX = strip_matrix_OneBump(NStrip, l, W0) 	# Create connection matrix for the Strip cell network
MY = strip_matrix_OneBump(NStrip, l, W0) 	# Create connection matrix for the Strip cell network
connect_with_matrix(StripX, StripX, MX) 	# Apply the connection matrix to the Stripcell network
connect_with_matrix(StripY, StripY, MY) 	# Apply the connection matrix to the Stripcell network

# Connect strip cells to grid cells
MG = gridcell_matrix_from_updrl(NStrip)
connect_with_matrix(Strip, Grid, MG)

# Set I parameter
for i in range(0,2*NStrip):
	StripX[i].I = I
	StripY[i].I = I




# INITIAL CONDITIONS .... SET INITIAL TO THE ORIGIN
Right[NStrip/2].sp = 1.0
Left[NStrip/2].sp  = 1.0
Up[NStrip/2].sp    = 1.0
Down[NStrip/2].sp  = 1.0


t= 0
m = 0
d = 0
# X = arange(0,Nstrip)
# XX, YY = meshgrid(X,X)

sp = get_spike_rates(Grid)
sp = sp.reshape([NStrip,NStrip])
#fig = mlab.figure(size=(800,800))
#MAVI = mlab.surf(sp)
#SMAVI = MAVI.mlab_source
plt.figure()
plt.ion()



# MAIN TIMELOOP
while(t < 4000):
	t= t+dt
	m= m+1
	stepNetwork(Strip, t, dt)      # Perform time-step in Strip cells
	stepNetwork(Grid, t, dt)      	# Perform time-step in Strip cells
	updateNetwork(Strip)	    # Update Neuron.sp = Neuron.s
	updateNetwork(Grid)	     	# Update Neuron.sp = Neuron.s
	
	
	if(m%10 == 0):
		plt.clf()

		# Plot strip cells
		ll = (get_spike_rates(Left), get_spike_rates(Down))
		rr = (get_spike_rates(Right), get_spike_rates(Up))
		
		plt.subplot(2,2,1)
		plt.plot(xStrip, ll[0], '-o')
		plt.plot(xStrip, rr[0], '-s')
		plt.ylim((0,1))
		plt.title('Strip cell module 1')
		plt.subplot(2,2,4)
		plt.plot(ll[1], xStrip, '-o')
		plt.plot(rr[1], xStrip, '-s')
		plt.xlim((0,1))
		plt.title('Strip cell module 2')
			
		# Plot grid cells
		sp = get_spike_rates(Grid)
		sp = sp.reshape([NStrip,NStrip])
		plt.subplot(2,2,3)
		plt.pcolor(xGrid, yGrid, sp)
		plt.title('Grid cells')
#		SMAVI.reset(scalars=sp)
#		mlab.savefig('./fig/%d.png'%d)

		# Plot representation of space
		plt.subplot(2,2,2)
		plt.plot(position[0], position[1], 'o')
		plt.xlim(posXLim)
		plt.ylim(posYLim)

		plt.draw()
		
		d= d+1
		
	# Update velocity in circular motion
#	xvel = cos(t/20.0)
#	yvel = sin(t/20.0)
#	Head[0].sp = fmax(xvel, 0.0)
#	Head[1].sp = fmax(-xvel, 0.0)
#	Head[2].sp = fmax(yvel, 0.0)
#	Head[3].sp = fmax(-yvel, 0.0)

	# Update velocity randomly
	# Compute change in velocity randomly
	dv = array((velocity[1], -velocity[0]))*(2*random.random() - 1)
	dv = dv/linalg.norm(dv)*accel*dt
	# Update velocity and normalize
	velocity = velocity + dv
	velocity = velocity / linalg.norm(velocity)
	# Update position
	position -= velocity*dt*0.1
	# print position
	# print("New velocity: (%f, %f)" % (velocity[0], velocity[1]))
	# Update HD cells
	Head[0].sp = fmax(velocity[0], 0.0)
	Head[1].sp = fmax(-velocity[0], 0.0)
	Head[2].sp = fmax(velocity[1], 0.0)
	Head[3].sp = fmax(-velocity[1], 0.0)
