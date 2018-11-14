#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


with open('totals.dat','r') as f:
    raw_times = f.readline().strip()
    times = np.array([float(time) for time in raw_times.split('\t')])
    raw_left = f.readline().strip()
    left = np.array([float(left) for left in raw_left.split('\t')])
    raw_right = f.readline().strip()
    right = np.array([float(right) for right in raw_right.split('\t')])
print(left)
plt.plot(times[::100],right[::100])
plt.plot(times[::100],left[::100])
plt.show()
plt.plot(times[::100],np.abs(left - right)[::100])
plt.show()
