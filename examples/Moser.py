import sys
sys.path.append('../')

from numpy import *
from nlab import *
import matplotlib.pyplot as plt
from mayavi import mlab

NGrid = 2**10
Nx = sqrt(NGrid)
dt = 0.1

alpha = 0.1
l  = 3
I = 3.0


# The phi function
def phi(x1, y1, x2, y2, u, v):
	l = 2
	R = 15
	W0 = -0.2
	dst = sqrt((x1-x2-l*u)**2 + (y1-y2-l*v)**2)
	
	if(dst < R):
		return W0
	else:
		return 0




# Create Head cells
Head = array([Neuron_IF() for _ in range(0,4)])

# Create Grid cells
Grid  = array([Neuron_IF() for _ in range(0,4*NGrid)])
Up    = Grid[0:NGrid]
Down  = Grid[NGrid:2*NGrid]
Right = Grid[2*NGrid:3*NGrid]
Left  = Grid[3*NGrid:4*NGrid]


# Set Input factor I
for i in range(0,4*NGrid):
	Grid[i].I = I

# Connect Head Cells
connect_one_to_many(Head[0], Up, alpha)
connect_one_to_many(Head[1], Down, alpha)
connect_one_to_many(Head[2], Right, alpha)
connect_one_to_many(Head[3], Left, alpha)


# Create connectivity
print 'Generating connectivity (might take a long time)'
M = gridcell_matrix_from_phi(NGrid, phi)
print 'Applying connectivity'
connect_with_matrix(Grid, Grid, M)


# Initialize
P = NGrid/2 + Nx/2
Up[P].sp = 1.0
Down[P].sp = 1.0
Left[P].sp = 1.0
Right[P].sp = 1.0



t= 0
m = 0
d = 0
X = arange(0,Nx)
XX, YY = meshgrid(X,X)

sp = get_spike_rates(Up)
sp = sp.reshape([Nx,Nx])
fig = mlab.figure(size=(800,800))
MAVI = mlab.surf(sp)
SMAVI = MAVI.mlab_source


print 'Starting main loop'
# MAIN TIMELOOP
while(t < 1000):
	t= t+dt
	m= m+1
	stepNetwork(Grid, dt)      	# Perform time-step in Strip cells
	updateNetwork(Grid)	     	# Update Neuron.sp = Neuron.s
	
	
	if(m%60 == 0):
		sp = get_spike_rates(Up)
		sp = sp.reshape([Nx,Nx])
		SMAVI.reset(scalars=sp)
		mlab.savefig('./fig/%d.png'%d)
		d= d+1
		
	# xvel = cos(t/20.0)
	# yvel = sin(t/20.0)
	# Head[0].sp = fmax(xvel, 0.0)
	# Head[1].sp = fmax(-xvel, 0.0)
	# Head[2].sp = fmax(yvel, 0.0)
	# Head[3].sp = fmax(-yvel, 0.0)
	if(m==80):
		Head[0].sp = 1.0
	if(m==1000):
		Head[1].sp = 1.0
		Head[0].sp = 0.0
		
	
