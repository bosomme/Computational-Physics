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
filename = "Solar_system.txt"
#Change this name!#

if filename == ".txt":
   sys.exit("Error: please put in filename, title and labels")

file = open(filename, 'r')
file.readline()

data = file.readlines()
x_sun = [];     y_sun = []
x_mercury = []; y_mercury = []
x_venus = [];   y_venus = []
x_earth = [];   y_earth = []
x_mars = [];    y_mars = []
x_jupiter = []; y_jupiter = []
x_saturn = [];  y_saturn = []
x_uranus = [];  y_uranus = []
x_neptune = []; y_neptune = []


for line in data:
    x_sun.append(line.split()[0]);      y_sun.append(line.split()[1])
    x_mercury.append(line.split()[2]);  y_mercury.append(line.split()[3])
    x_venus.append(line.split()[4]);    y_venus.append(line.split()[5])
    x_earth.append(line.split()[6]);    y_earth.append(line.split()[7])
    x_mars.append(line.split()[8]);     y_mars.append(line.split()[9])
    x_jupiter.append(line.split()[10]); y_jupiter.append(line.split()[11])
    x_saturn.append(line.split()[12]);  y_saturn.append(line.split()[13])
    x_uranus.append(line.split()[14]);  y_uranus.append(line.split()[15])
    x_neptune.append(line.split()[16]); y_neptune.append(line.split()[17])


file.close()

fig = plt.figure(figsize=(6,4.5))

ax1 = fig.add_subplot(111)

ax1.set_xlabel('x-value [AU]')
ax1.set_ylabel('y-value [AU]')

ax1.plot(x_sun, y_sun)
ax1.plot(x_mercury, y_mercury)
ax1.plot(x_venus, y_venus)
ax1.plot(x_earth, y_earth)
ax1.plot(x_mars, y_mars)
ax1.plot(x_jupiter, y_jupiter)
ax1.plot(x_saturn, y_saturn)
ax1.plot(x_uranus, y_uranus)
ax1.plot(x_neptune, y_neptune)
#ax1.plot(x[-1], y[-1], 'rx')
#ax1.annotate('Starting Point', xy=(0.5,0), textcoords='data')
ax1.axis('equal')   #Equal axis

leg = ax1.legend(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'])

plt.show()
