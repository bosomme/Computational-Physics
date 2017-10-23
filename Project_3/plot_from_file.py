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

#Change this name!#
filename = "escape_velocity.txt"
#Change this name!#

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



fig = plt.figure(figsize=(6,4.5))

ax1 = fig.add_subplot(111)

ax1.set_xlabel('x-value [AU]')
ax1.set_ylabel('y-value [AU]')

ax1.plot(x, y)
ax1.plot(x[-1], y[-1], 'rx')
#ax1.annotate('Starting Point', xy=(0.5,0), textcoords='data')
ax1.axis('equal')   #Equal axis

leg = ax1.legend(['Earth'])

plt.show()
