#!/usr/bin/env python3
import numpy as np
import json
import sys


def double_channel_param_set(param_coupling):
    resolution = 25
    for option, option_values in param_coupling['options'].items():
        num_particles = option_values['num_particles']
        identifier = "double_channel_data/params_%s" % num_particles
        param_coupling['defaults']['num_particles'] = num_particles
        with open(identifier + '.in', 'w') as f:
            for initial_ratio in param_coupling['initial_ratio']:
                param_coupling['defaults']['initial_ratio'] = initial_ratio
                for name, interval in param_coupling['ranges'].items():
                    values = param_coupling['defaults'].copy()
                    values['threshold'] *= option_values['threshold_multiplier']
                    sampled_range = np.linspace(interval[0], interval[1], resolution)
                    for i in range(resolution):
                        values[name] = sampled_range[i]
                        if name == 'threshold':
                            # Bit of a dirty hack because I don't want to
                            # make the structure of the parameter dependent on the number of particles
                            values[name] *= option_values['threshold_multiplier']
                        cmd = ["%.4f" % values[p_name] for p_name in values.keys()] + [identifier]
                        f.write(" ".join(cmd) + "\n")


def double_channel_heat_map(param_coupling):
    resolution = 50
    for option, option_values in param_coupling['options'].items():
        num_particles = option_values['num_particles']
        identifier = "double_channel_data/heatmap_%s" % num_particles
        param_coupling['defaults']['num_particles'] = num_particles
        with open(identifier + '.in', 'w') as f:
            t_bounds = param_coupling['ranges']['threshold']
            threshold_interval = np.linspace(t_bounds[0], t_bounds[1], resolution) * option_values[
                'threshold_multiplier']
            w_bounds = param_coupling['ranges']['second_width']
            width_interval = np.linspace(w_bounds[0], w_bounds[1], resolution)
            for initial_ratio in [0.25, 0.5]:
                param_coupling['defaults']['initial_ratio'] = initial_ratio
                values = param_coupling['defaults'].copy()
                for i in range(resolution):
                    for j in range(resolution):
                        values['threshold'] = threshold_interval[i]
                        values['second_width'] = width_interval[j]
                        cmd = ["%.4f" % values[p_name] for p_name in values.keys()] + [identifier]
                        f.write(" ".join(cmd) + "\n")


if __name__ == '__main__':
    filename = 'params_double_channel.json'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    with open(filename, 'r') as param_file:
        parameter_sets = json.load(param_file)
        double_channel_param_set(parameter_sets)
        double_channel_heat_map(parameter_sets)
