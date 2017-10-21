# -*- coding: utf-8 -*-
"""
Program that plots from .txt-file
First column is x-value
Second column is y-value

Put name of txt-file in "filename", and edit title and labels

Always reset filename = ".txt"
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = "velocity_verlet.txt"

if filename == ".txt":
   sys.exit("Error: please put in filename, title and labels")

file = open(filename, 'r')
file.readline()

data = file.readlines()
t = []; x = []; y = []

for line in data:
    t.append(line.split()[0])
    x.append(line.split()[1])
    y.append(line.split()[2])

file.close()

x_init = x[0]; y_init = y[0]



fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_xlabel('x-value [AU]')
ax1.set_ylabel('y-value [AU]')

ax1.plot(x, y)
ax1.annotate('Starting Point', xy=(x[0],y[0]), textcoords='data')

leg = ax1.legend(['Earth, Orbit'])

plt.show()
