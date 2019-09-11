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
                cmd = ["%.2f" % values[name] for name in default.keys()] + [identifier]
                cmd_string = " ".join(cmd)
                f.write(cmd_string + "\n")

parameter_sets = \
    {"small": {"relations": [{"channel_length": [0.2, 0.8], "channel_width": [0.3, 1]},
                             {"channel_width": [0.4, 1.4], "urn_radius": [0.8, 1.5]},
                             {"channel_length": [0.1, 1], "urn_radius": [0.7, 1.3]},
                             {"channel_width": [0.3, 1], "threshold": [4, 22]},
                             {"channel_length": [0.3, 1], "threshold": [4, 22]},
                             {"urn_radius": [0.7, 1.3], "threshold": [4, 22]}],
               "defaults": {"channel_length": 0.5, "channel_width": 0.5, "urn_radius": 1, "threshold": 10,
                            "num_particles": 1000, "M_t": 1E5, "M_f": 1.5E5}},
     "medium": {"relations": [{'channel_length': [0.2, 0.8], 'channel_width': [0.3, 1]},
                              {'channel_width': [0.4, 1.4], 'urn_radius': [0.8, 1.5]},
                              {'channel_length': [0.1, 1], 'urn_radius': [0.7, 1.3]},
                              {'channel_width': [0.3, 1], 'threshold': [10, 64]},
                              {'channel_length': [0.3, 1], 'threshold': [10, 64]},
                              {'urn_radius': [0.7, 1.3], 'threshold': [10, 64]}],
                "defaults": {"channel_length": 0.5, "channel_width": 0.5, "urn_radius": 1, "threshold": 30,
                             "num_particles": 3000, "M_t": 5E5, "M_f": 6E5}},
     "large": {"relations": [{'channel_length': [0.2, 0.8], 'channel_width': [0.3, 1]},
                             {'channel_width': [0.4, 1.4], 'urn_radius': [0.8, 1.5]},
                             {'channel_length': [0.1, 1], 'urn_radius': [0.7, 1.3]},
                             {'channel_width': [0.3, 1], 'threshold': [40, 220]},
                             {'channel_length': [0.3, 1], 'threshold': [40, 220]},
                             {'urn_radius': [0.3, 1], 'threshold': [40, 220]}],
               "defaults": {"channel_length": 0.5, "channel_width": 0.5, "urn_radius": 1, "threshold": 300,
                            "num_particles": 10000, "M_t": 1E5, "M_f": 1.5E5}}}

if __name__ == '__main__':
    for size in parameter_sets.keys():
        default_values = parameter_sets[size]['defaults']
        for i, param in enumerate(parameter_sets[size]['relations']):
            generate_param_set(param, default_values, "param_file_%s_%d" % (size,i))
