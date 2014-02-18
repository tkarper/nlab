import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt
from mayavi import mlab

NStrip = 2**5
dt = 0.1

alpha = 0.1
l  = 3
W0 = -0.5
I = 1.0

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
fig = mlab.figure(size=(800,800))
MAVI = mlab.surf(sp)
SMAVI = MAVI.mlab_source



# MAIN TIMELOOP
while(t < 4000):
	t= t+dt
	m= m+1
	stepNetwork(Strip, t, dt)      # Perform time-step in Strip cells
	stepNetwork(Grid, t, dt)      	# Perform time-step in Strip cells
	updateNetwork(Strip)	    # Update Neuron.sp = Neuron.s
	updateNetwork(Grid)	     	# Update Neuron.sp = Neuron.s
	
	
	if(m%60 == 0):
		sp = get_spike_rates(Grid)
		sp = sp.reshape([NStrip,NStrip])
		SMAVI.reset(scalars=sp)
		mlab.savefig('./fig/%d.png'%d)
		d= d+1
		
	xvel = cos(t/20.0)
	yvel = sin(t/20.0)
	Head[0].sp = fmax(xvel, 0.0)
	Head[1].sp = fmax(-xvel, 0.0)
	Head[2].sp = fmax(yvel, 0.0)
	Head[3].sp = fmax(-yvel, 0.0)
	#if(m==80):
	#	Head[0].sp = 1.0
	#if(m==1000):
	#	Head[1].sp = 1.0
	#	Head[0].sp = 0.0
		
	
