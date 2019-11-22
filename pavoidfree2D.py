import numpy as np
from numpy.linalg import inv
import pylab, random
import matplotlib.pyplot as plt

from math import factorial

T  = 1024 # T should be even
Deltax = 4 # Deltax should be even
trials = 1000000
navoid = 0
for i in range(trials):
	poly1 = [ 1 for i in range(T/2) ]
	polym1 = [ -1 for i in range(T/2) ]
	poly = poly1 + polym1
	Polya = np.random.permutation(poly)
	Polyb = np.random.permutation(poly)
	#building the paths and verifying if they avoid
	Patha = [0 for i in range(T+1)]
	Pathb = [0 for i in range(T+1)]
	Patha[0] = -Deltax//2
	Pathb[0] = Deltax//2
	time = 1
	Patha[time] = Patha[time-1] + Polya[time-1]
	Pathb[time] = Pathb[time-1] + Polyb[time-1]
	while (Patha[time] != Pathb[time] and time != T):
		time = time+1
		Patha[time] = Patha[time-1] + Polya[time-1]
		Pathb[time] = Pathb[time-1] + Polyb[time-1]
	if time == T:
		navoid = navoid + 1
probavoid = float(navoid)/trials
print("T= ", T, " Deltax= ", Deltax, " probavoid = ", probavoid)
	
		
