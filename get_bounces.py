#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv('bounces.dat',header=None,names=['x','y'],sep=' ',nrows=5000)
df.plot.scatter('x','y')
plt.show()
