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
filename = "Earth_Jupiter.txt"
#Change this name!#

if filename == ".txt":
   sys.exit("Error: please put in filename, title and labels")

file = open(filename, 'r')
file.readline()

data = file.readlines()
x_earth = []; y_earth = []
x_jupiter = []; y_jupiter = []
x_sun = []; y_sun = []

for line in data:
    x_earth.append(line.split()[0])
    y_earth.append(line.split()[1])
    x_jupiter.append(line.split()[2])
    y_jupiter.append(line.split()[3])
    x_sun.append(line.split()[4])
    x_sun.append(line.split()[5])

file.close()



fig = plt.figure(figsize=(6,4.5))

ax1 = fig.add_subplot(111)

ax1.set_xlabel('x-value [AU]')
ax1.set_ylabel('y-value [AU]')

ax1.plot(x_earth, y_earth)
ax1.plot(x_jupiter, y_jupiter)
#ax1.plot(x_sun, y_sun)
#ax1.plot(x[-1], y[-1], 'rx')
#ax1.annotate('Starting Point', xy=(0.5,0), textcoords='data')
ax1.axis('equal')   #Equal axis

leg = ax1.legend(['Earth', 'Jupiter'])

plt.show()
