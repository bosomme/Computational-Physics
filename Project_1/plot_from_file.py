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

filename = ".txt"

if filename == ".txt":
   sys.exit("Error: please put in filename, title and labels")

file = open(filename, 'r')
file.readline()

data = file.readlines()
x = []; y1 = []; y2 = []

for line in data:
    x.append(line.split()[0])
    y1.append(line.split()[1])
    y2.append(line.split()[2])

file.close()


fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title('Special 1000x1000')
ax1.set_xlabel('x-value')
ax1.set_ylabel('y-value')

ax1.plot(x,y1,x,y2)

leg = ax1.legend(['Numerical','Analytical'])

plt.show()
