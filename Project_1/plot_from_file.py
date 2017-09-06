# -*- coding: utf-8 -*-
"""
Program that plots from .txt-file
First column is x-value
Second column is y-value

Put name of txt-file in "filename", and edit title and labels
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = "test1.txt"

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

ax1.set_title('Plot title...')
ax1.set_xlabel('your x label..')
ax1.set_ylabel('your y label...')

ax1.plot(x,y1,x,y2)

leg = ax1.legend(['Numerical','Analytical'])

plt.show()
