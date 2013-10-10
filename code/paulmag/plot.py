import matplotlib.pyplot as plt
import numpy as np

infile = open("data/prob_no_coulomb.dat", "r") # No coulomb interaction.
#infile = open("data/prob_coulomb0.dat", "r") # Coulomb interaction with five
#infile = open("data/prob_coulomb1.dat", "r") # different oscillator strengths.
#infile = open("data/prob_coulomb2.dat", "r")
#infile = open("data/prob_coulomb3.dat", "r")
#infile = open("data/prob_coulomb4.dat", "r")

info = infile.readline().split(",")
n, rhoMax, omega_r = int(info[0]), float(info[1]), float(info[2])
rho = np.linspace(0, rhoMax, n);

wave = [ [], [], [] ]

for  line in infile:
    line = line.split(",")
    for k in range(3):
        wave[k].append(float( line[k] ))


wave = np.array(wave)

for k in range(3):
    wave[k, 1:] /= rho[1:]
    wave[k, 0]  /= rho[1] * 0.5 # avoid zero division (not neccesarily accurate)

for k in range(3):
    for i in range(n):
        wave[k,i] = wave[k,i]**2 # find probability
    wave[k] = n * wave[k] / np.sum(wave[k]) / rhoMax # normalize

plt.plot(rho,wave[0], rho,wave[1], rho,wave[2])

plt.xlabel("rho")
plt.ylabel("probability")

plt.title("Two electrons in oscillator, no coulomb force")
#plt.title("Two electrons in oscillator, coulomb force, omega_r=%g" % omega_r)

plt.legend(("$\lambda_0$", "$\lambda_1$", "$\lambda_2$"), loc="best")

plt.show()
