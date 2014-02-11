import sys
sys.path.append('../')


from numpy import *
from nlab import *
import matplotlib.pyplot as plt

n = 40
L = array([Neuron_IF() for _ in range(0,n)])
#for i in range(0,n):
#	L[i].sp = i
	
connect_one_to_many(L[0], L, 1.0)
