#!/usr/bin/env python3
import numpy as np
import json
import sys


def double_channel_param_set(param_coupling):
    resolution = 5
    for num_particles in param_coupling['options']['num_particles']:
        identifier = "double_channel_data/params_%s" % num_particles
        param_coupling['defaults']['num_particles'] = num_particles
        with open(identifier + '.in', 'w') as f:
            for initial_ratio in param_coupling['options']['initial_ratio']:
                param_coupling['defaults']['initial_ratio'] = initial_ratio
                for name, interval in param_coupling['ranges'].items():
                    values = param_coupling['defaults'].copy()
                    sampled_range = np.linspace(interval[0], interval[1], resolution)
                    for i in range(resolution):
                        values[name] = sampled_range[i]
                        cmd = ["%.3f" % values[pname] for pname in values.keys()] + [identifier]
                        f.write(" ".join(cmd) + "\n")


if __name__ == '__main__':
    filename = 'params_double_channel.json'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    with open(filename, 'r') as param_file:
        parameter_sets = json.load(param_file)
        double_channel_param_set(parameter_sets)
