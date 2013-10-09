import matplotlib.pyplot as plt
import numpy as np

"""
infile0 = open("data/prob_no_coulomb_0.dat", "r")
infile1 = open("data/prob_no_coulomb_1.dat", "r")
infile2 = open("data/prob_no_coulomb_2.dat", "r")
"""
infile0 = open("data/prob_coulomb2_0.dat", "r")
infile1 = open("data/prob_coulomb2_1.dat", "r")
infile2 = open("data/prob_coulomb2_2.dat", "r")


wave0 = []
wave1 = []
wave2 = []

n = 0

for line0,line1,line2 in zip(infile0, infile1, infile2):
    wave0.append(float(line0))
    wave1.append(float(line1))
    wave2.append(float(line2))
    n += 1

r = np.linspace(0, 1, n);

# Normalize probabilities:
wave0, wave1, wave2 = np.array(wave0), np.array(wave1), np.array(wave2)
wave0 = n * wave0 / np.sum(wave0)
wave1 = n * wave1 / np.sum(wave1)
wave2 = n * wave2 / np.sum(wave2)

plt.plot(r,wave0, r,wave1, r,wave2)
plt.xlabel("r")
plt.ylabel("probability")
plt.show()
