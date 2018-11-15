#!/usr/bin/env python3
import numpy as np
import subprocess

import scipy.io as sio

radii = np.linspace(0.1,1.5,15)
capacities = np.arange(1,19,2)
results = np.zeros([len(radii),len(capacities)])
for i in range(len(radii)):
    for j in range(len(capacities)):
        radius = radii[i]
        capacity = capacities[j]
        popen = subprocess.check_output(['./terrier',str(radius),str(capacity)])
        results[i,j] = float(popen)
        print("For radius %.2f and capacity %d we have thermalization time %.2f"%(radius,capacity,results[i,j]))
sio.savemat('therm',{'radii':radii,'capacities':capacities,'results':results})