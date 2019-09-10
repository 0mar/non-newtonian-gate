#!/usr/bin/env python3
import numpy as np


def generate_param_set(param_coupling, default, identifier='param_file1'):
    values = default.copy()
    res = 10
    it = iter(param_coupling)
    x = next(it)
    x_vals = np.linspace(*param_coupling[x], res)
    y = next(it)
    y_vals = np.linspace(*param_coupling[y], res)
    with open(identifier + ".in", 'w') as f:
        # f.write("# Parameters: %s \n" % " | ".join(param_names))
        for x_val in x_vals:
            values[x] = x_val
            for y_val in y_vals:
                values[y] = y_val
                cmd = ["%.2f" % values[param] for param in param_names] + [identifier]
                cmd_string = " ".join(cmd)
                f.write(cmd_string + "\n")


param_names = ['channel_width', 'channel_length', 'urn_radius', 'threshold']
params = [{'channel_length': [0.2, 0.8], 'channel_width': [0.3, 1]},
          {'channel_width': [0.4, 1.4], 'urn_radius': [0.8, 1.5]},
          {'channel_length': [0.1, 1], 'urn_radius': [0.7, 1.3]},
          {'channel_width': [0.3, 1], 'threshold': [10, 65]},
          {'channel_length': [0.3, 1], 'threshold': [10, 65]},
          {'urn_radius': [0.3, 1], 'threshold': [10, 65]}]
default_values = {'channel_width': 0.5, 'channel_length': 0.5, 'urn_radius': 1, 'threshold': 10}

if __name__ == '__main__':
    for i, param in enumerate(params):
        generate_param_set(param, default_values, "param_file_%d" % i)
