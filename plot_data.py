#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

plot_dir = 'plots'
single_channel_dir = 'single_channel_data'
double_channel_dir = 'double_channel_data'


def threshold_function(params):
    length, width, radius, threshold, num_particles = params['length'], params['width'], params['radius'], params[
        'threshold'], params['num_particles']
    area = np.pi * radius ** 2 - radius ** 2 * np.arcsin(width / (2 * radius)) + width * np.sqrt(
        4 * radius ** 2 - width ** 2) / 4
    return num_particles * width * length / (8 * threshold * area * 0.75)


def plot_single_channel_heat_map(filename, num_particles):
    param_names = ["length", "width", "radius", "threshold"]
    plt.figure(figsize=(22, 15))
    for i in range(6):
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

        ax = plt.subplot(2, 3, i + 1)
        plt.scatter(x, y, c=chi, s=150, alpha=0.7, cmap='PiYG')
        plt.colorbar()
        plt.contour(X, Y, Z, [1], linewidths=3, linestyles='dashdot', colors='black')

        plt.xlabel(x_label.title())
        plt.ylabel(y_label.title())
        if y_label == 'threshold':
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.savefig("%s/%s/omar-recreation-%s.pdf" % (single_channel_dir, plot_dir, num_particles))
    plt.close()


def plot_single_channel_heat_maps():
    for size in ["small", "large"]:
        filename = '%s/param_file_%s_%s.out' % (single_channel_dir, size, "%d")
        num_particles = {"small": "1E3", "large": "1E4"}[size]
        plot_single_channel_heat_map(filename, num_particles)


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
        # Todo: unique the data
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


def plot_double_channel_data_manuscript():
    df = pd.DataFrame()
    for num_particles in [1000, 10000]:
        file_id = '%s/params_%d' % (double_channel_dir, num_particles)
        try:
            part_df = pd.read_csv(file_id + '.out', header=None, sep=',',
                                  names=['Relative threshold', 'Width of second channel', 'Length of second channel',
                                         'initial_ratio', 'Mass spread', 'Relative current'])
            df = df.append(part_df.assign(num_particles=num_particles))
        except FileNotFoundError:
            print("%s not found, continuing" % file_id)
            continue

    measurements = ['Mass spread', 'Relative current']
    variables = {'Relative threshold', 'Width of second channel', 'Length of second channel'}
    markers = {0.25: 'v', 0.5: 's', 0.75: 'o'}
    marker_styles = {1000: 'none', 10000: 'full'}
    np_format = {1000: '10^3', 10000: '10^4'}
    df.loc[:, 'Relative threshold'] /= df.num_particles
    df.loc[:, 'Relative current'] /= df.num_particles
    for measurement in measurements:
        for variable in variables:
            legend = []
            plt.figure(figsize=(7, 7))
            other_variables = variables.copy()
            other_variables.remove(variable)
            for num_particles in [1000, 10000]:
                for initial_ratio in [0.25, 0.75]:
                    sdf = df[(df.initial_ratio == initial_ratio) & (df.num_particles == num_particles)].copy()
                    for other_variable in other_variables:
                        sdf = sdf[sdf[other_variable] == sdf[other_variable].value_counts().index[0]].sort_values(
                            variable)
                    plt.plot(sdf[variable].values, sdf[measurement].values, marker=markers[initial_ratio],
                             fillstyle=marker_styles[num_particles], color='black', linestyle='-.')
                    legend.append('$N=%s$ , $\\chi_0 = %.1f$' % (np_format[num_particles], initial_ratio * 2 - 1))
            plt.xlabel(variable)
            plt.ylabel(measurement)
            plt.legend(legend)
            plt.savefig("%s/%s/%s_vs_%s.pdf" % (
                double_channel_dir, plot_dir, variable.lower().replace(' ', '_'),
                measurement.lower().replace(' ', '_')))
            plt.close()


def plot_double_channel_heatmap():
    for num_particles in [1000, 10000]:
        filename = '%s/heatmap_%d.out' % ('double_channel_data', num_particles)
        outputs = ['chi', 'current']
        df = pd.read_csv(filename, header=None,
                         names=["threshold", "second_width", "second_length", "initial_ratio", "chi", "current"])
        df.loc[:, 'current'] = np.abs(df.current / num_particles)
        df.loc[:, 'chi'] = np.abs(df.chi)
        plt.figure(figsize=(11, 10))
        counter = 0
        for output in outputs:
            for init_ratio in df.initial_ratio.unique():
                sdf = df[df.initial_ratio == init_ratio]
                counter += 1
                plt.subplot(2, 2, counter)
                plt.scatter(sdf.threshold.values / num_particles, sdf.second_width.values, c=sdf[output].values, s=100,
                            marker='o')
                plt.axis([0, 20 / 1000, 0, 0.1])
                plt.xlabel("Relative threshold")
                plt.ylabel("Second channel width")
                plt.title("%s, initial ratio = %.2f" % (output.title(), init_ratio))
                plt.colorbar()
        plt.savefig("double_channel_data/plots/double_channel_heatmap_%d.pdf" % num_particles)
        plt.show()


# plot_double_channel_data()
# plot_double_channel_data_manuscript()
plot_single_channel_heat_maps()
# plot_double_channel_heatmap()
