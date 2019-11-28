#!/usr/bin/env python3
import numpy as np
import json
import sys

def single_channel_param_set(param_coupling, default, identifier):
    values = default.copy()
    res = 20
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

if __name__ == '__main__':
    filename = 'params_single_channel.json'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    with open(filename,'r') as param_file:
        parameter_sets = json.load(param_file)
    for size in parameter_sets.keys():
        default_values = parameter_sets[size]['defaults']
        for i, param in enumerate(parameter_sets[size]['relations']):
            single_channel_param_set(param, default_values, "single_channel_data/param_file_%s_%d" % (size,i))
