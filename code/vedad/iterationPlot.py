import matplotlib.pyplot as plt
import numpy as np

def read():

	inFile = open("data/iterationsData.txt",'r')

	global _n; _n = []
	global _iterations; _iterations = []

	for line in inFile:
		columns = line.split()
		_n.append(float(columns[0]))
		_iterations.append(float(columns[1]))

	inFile.close()

read()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(_n,_iterations,'o', linestyle='-')
fig.suptitle('Number of iterations as function of\n dimensionality of matrix')
ax.set_xlabel('$n$', fontsize='16')
ax.set_ylabel('Iterations')
plt.grid('on')
plt.show()
