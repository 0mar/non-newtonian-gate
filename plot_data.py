#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

plot_dir = 'plots'
single_channel_dir = 'single_channel_data'


def threshold(params):
    length, width, radius, threshold, num_particles = params['length'], params['width'], params['radius'], params[
        'threshold'], params['num_particles']
    area = np.pi * radius ** 2 - radius ** 2 * np.arcsin(width / (2 * radius)) + width * np.sqrt(
        4 * radius ** 2 - width ** 2) / 4
    return num_particles*width*length/(8*threshold*area*0.75)


def plot(filename, num_particles):
    param_names = ["length","width","radius","threshold"]
    plt.figure(figsize=(22,15))
    for i in range(6):
        row = i//3
        col = i%3
        df = pd.read_csv(filename % i, header=None, names=["length", "width", "radius", "threshold", "chi"])
        df['num_particles']=int(float(num_particles))
        sub_df = df.loc[:, (df != df.iloc[0]).any()]
        x_label = sub_df.columns[0]
        y_label = sub_df.columns[1]
        x = df[x_label].values
        y = df[y_label].values
        chi = df[sub_df.columns[2]].values
        unused_cols = set(param_names) - set(sub_df.columns)
        
        param1_name = unused_cols.pop()
        param2_name = unused_cols.pop()
        param1_val = df[param1_name][0]
        param2_val = df[param2_name][0]
        params = ",".join(param_names+["num_particles"]).replace(x_label,'x').replace(y_label,'y')
       
        x_ = np.linspace(np.min(x),np.max(x))
        y_ = np.linspace(np.min(y),np.max(y))
        X,Y = np.meshgrid(x_,y_)
        Z = threshold({param1_name:param1_val,param2_name:param2_val,x_label:X,y_label:Y,'num_particles':int(float(num_particles))})
        
        plt.subplot(2,3,i+1)
        plt.scatter(x,y,c=chi,s=300,marker='s',cmap='autumn')
        plt.colorbar()
        plt.contour(X,Y,Z,[1])
        
        plt.xlabel(x_label)
        plt.ylabel(y_label)
    plt.savefig("%s/%s/omar-recreation-%s.pdf" % (single_channel_dir, plot_dir, num_particles))

for size in ["small","large"]:
    filename = '%s/param_file_%s_%s.out'%(single_channel_dir,size,"%d")
    num_particles = {"small":"1E3","large":"1E4"}[size]
    plot(filename, num_particles)

