#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

plot_dir = 'input/plots'
input_dir = 'input'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

def plot(filename,size):
    plt.figure(figsize=(22,15))
    for i in range(6):
        row = i//3
        col = i%3
        df = pd.read_csv(filename%i,header=None,names=["length","width","radius","threshold","chi"])
        df = df.loc[:, (df!= df.iloc[0]).any()]
        x = df[df.columns[0]].values
        y = df[df.columns[1]].values
        chi = df[df.columns[2]].values
        plt.subplot(2,3,i+1)
        plt.scatter(x,y,c=chi,s=500)
        plt.colorbar()
        plt.xlabel(df.columns[0])
        plt.ylabel(df.columns[1])
    plt.savefig("%s/omar-recreation-%s.pdf"%(plot_dir,size))

for size in ["small","medium","large"]:
    filename = '%s/param_file_%s_%s.out'%(input_dir,size,"%d")
    plot(filename,size)

intro = """
reset
set notitle
set nokey
set size square

set view map

set contour
set cntrparam level incremental 1., 1., 1.
#unset contour
#set parametric
set dgrid3d
set pm3d interpolate 0,0
set hidden3d
set xlabel font ",17"
set ylabel font ",17"
set ylabel offset -3,0
set tics font ",13"
num_particles={num_particles:.2f}
area(L,W,R,N)=pi*R**2-R**2*asin(W/(2*R))+W*sqrt(4*R**2-W**2)/4.

threshold_line(L,W,R,T,N)=N*W*L/(8.*T*area(L,W,R,N))

"""

plot_string = """
#Graph {x_label}-{y_label}
set xrange [{x1}:{x2}]
set yrange [{y1}:{y2}]
set cbrange [0:1]
set xlabel "{x_label}"
set ylabel "{y_label}"

{param1_name}={param1_value:.2f}
{param2_name}={param2_value:.2f}
f(x,y)=threshold_line({params})
set terminal pdf
set output "{plot_name}"
splot '{out_name}' u {ratio} with pm3d nocontour, f(x,y) with lines lw 3 dt 10 lc rgb "white"

"""

def output_gnu(prefix=''):
    filename_format='param_file_%s_%d.out'
    nums_particles=["1E3","3E3","1E4"]
    data_name_fmt='%s/dyn-N%s-%s-%s.dat'
    param_names = ["length","width","radius","threshold"]
    ratios = ['4:5:1','5:3:1','4:3:1','5:6:1','4:6:1','3:6:1']
    for si,size in enumerate(["small","medium","large"]):
        num_particles = nums_particles[si]
        gnu_plot_filename = '%s/gnu_plot_%s.gpi'%(input_dir,size)
        with open(gnu_plot_filename,'w') as gnu_plot_file:
            gnu_plot_file.write(intro.format(num_particles=int(float(num_particles))))
            for i in range(6):
                df = pd.read_csv(filename_format%(size,i),header=None,names=["length","width","radius","threshold","chi"])
                df['num_particles']=int(float(num_particles))
                sub_df = df.loc[:, (df != df.iloc[0]).any()]
                x_label = sub_df.columns[0]
                y_label = sub_df.columns[1]
                data_name = data_name_fmt%(input_dir/num_particles,x_label,y_label)
                df.to_csv(data_name,sep='\t',columns=['chi','num_particles','radius','length','width','threshold'], index=False, header=False)
                unused_cols = set(param_names) - set(sub_df.columns)
                x1 = df[x_label].min()
                x2 = df[x_label].max()
                y1 = df[y_label].min()
                y2 = df[y_label].max()
                param1_name = unused_cols.pop()
                param2_name = unused_cols.pop()
                param1_val = df[param1_name][0]
                param2_val = df[param2_name][0]
                params = ",".join(param_names+["num_particles"]).replace(x_label,'x').replace(y_label,'y')

                plot_name = "%s/plot-N%s-%s-%s.pdf"%(plot_dir,num_particles, x_label,y_label)
                full_string = plot_string.format(x_label=x_label, y_label=y_label,x1=str(x1),x2=str(x2),y1=str(y1),y2=str(y2),
                                                 param1_name=param1_name,param2_name=param2_name,param1_value=param1_val,param2_value=param2_val,
                                                 params=params,plot_name=plot_name,out_name=data_name,ratio=ratios[i])
                gnu_plot_file.write(full_string)


output_gnu()



