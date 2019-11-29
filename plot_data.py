#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plot_dir = 'plots'
single_channel_dir = 'single_channel_data'
double_channel_dir = 'double_channel_data'


def threshold_function(params):
    length, width, radius, threshold, num_particles = params['length'], params['width'], params['radius'], params[
        'threshold'], params['num_particles']
    area = np.pi * radius ** 2 - radius ** 2 * np.arcsin(width / (2 * radius)) + width * np.sqrt(
        4 * radius ** 2 - width ** 2) / 4
    return num_particles * width * length / (8 * threshold * area * 0.75)


def plot(filename, num_particles):
    param_names = ["length", "width", "radius", "threshold"]
    plt.figure(figsize=(22, 15))
    for i in range(6):
        row = i // 3
        col = i % 3
        df = pd.read_csv(filename % i, header=None, names=["length", "width", "radius", "threshold", "chi"])
        df['num_particles'] = int(float(num_particles))
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
        params = ",".join(param_names + ["num_particles"]).replace(x_label, 'x').replace(y_label, 'y')

        x_ = np.linspace(np.min(x), np.max(x))
        y_ = np.linspace(np.min(y), np.max(y))
        X, Y = np.meshgrid(x_, y_)
        Z = threshold_function({param1_name: param1_val, param2_name: param2_val, x_label: X, y_label: Y,
                                'num_particles': int(float(num_particles))})

        plt.subplot(2, 3, i + 1)
        plt.scatter(x, y, c=chi, s=150, alpha =0.7, cmap='PiYG')
        plt.colorbar()
        plt.contour(X, Y, Z, [1], linewidths=3, linestyles='dashdot', colors='black')

        plt.xlabel(x_label)
        plt.ylabel(y_label)
    plt.savefig("%s/%s/omar-recreation-%s.pdf" % (single_channel_dir, plot_dir, num_particles))
    plt.close()


def plot_single_channel_data():
    for size in ["small", "large"]:
        filename = '%s/param_file_%s_%s.out' % (single_channel_dir, size, "%d")
        num_particles = {"small": "1E3", "large": "1E4"}[size]
        plot(filename, num_particles)


def plot_double_channel_data():
    for num_particles in [1000, 10000]:
        file_id = '%s/params_%d' % (double_channel_dir, num_particles)
        try:
            df = pd.read_csv(file_id + '.out', header=None, sep=',',
                             names=['threshold', 'second_width', 'second_length', 'initial_ratio', 'mass_spread',
                                    'current'])
        except FileNotFoundError:
            print("%s not found, continuing" % file_id)
            continue
        plt.figure(figsize=(5, 5))
        measurements = ['mass_spread', 'current']
        variables = {'threshold', 'second_width', 'second_length'}

        for measurement in measurements:
            for variable in variables:
                plt.figure()
                other_variables = variables.copy()
                other_variables.remove(variable)
                for initial_ratio in df.initial_ratio.unique():
                    sdf = df[df.initial_ratio == initial_ratio].copy()
                    for other_variable in other_variables:
                        sdf = sdf[sdf[other_variable] == sdf[other_variable].value_counts().index[0]].sort_values(
                            variable)
                    plt.plot(sdf[variable].values, sdf[measurement].values, 's-')
                    plt.xlabel(variable)
                plt.ylabel(measurement)
                plt.legend(df.initial_ratio.unique())
                plt.savefig(
                    "%s/%s/%s_to_%s_(%d).png" % (double_channel_dir, plot_dir, variable, measurement, num_particles))
                plt.close()


plot_double_channel_data()
plot_single_channel_data()
