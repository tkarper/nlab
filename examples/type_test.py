import sys
sys.path.append('../')


from numpy import *
from nlab import *
import matplotlib.pyplot as plt


L = array([Neuron_IF() for _ in range(0,4)])
for i in range(0,4):
	L[i].sp = i
	
connect_one_to_many(L[0], L, 1.0)
