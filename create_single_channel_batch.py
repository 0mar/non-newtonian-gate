#!/usr/bin/env python3
import numpy as np
import json
import sys


def single_channel_param_set(param_coupling, default, identifier):
    values = default.copy()
    res = 50
    it = iter(param_coupling)
    x = next(it)
    x_values = np.linspace(*param_coupling[x], res)
    y = next(it)
    y_values = np.linspace(*param_coupling[y], res)
    with open(identifier + ".in", 'w') as f:
        # f.write("# Parameters: %s \n" % " | ".join(param_names))
        for x_val in x_values:
            values[x] = x_val
            for y_val in y_values:
                values[y] = y_val
                cmd = ["%.4f" % values[name] for name in default.keys()] + [identifier]
                cmd_string = " ".join(cmd)
                f.write(cmd_string + "\n")


def single_channel_exploration_set(parameter_sets):
    resolution = 10
    relation_set = {}
    for size, value in parameter_sets.items():
        identifier = "single_channel_data/explorer_%s" % size
        for relation in value['relations']:
            relation_set.update(relation)
        relation_set.pop('threshold')
        intervals = []
        for relation in relation_set:
            intervals.append(np.linspace(relation_set[relation][0], relation_set[relation][1], resolution))
        if len(intervals) > 3:
            raise ValueError("Don't think this will work in higher than 3 dimensions atm")
        meshes = np.meshgrid(*intervals)
        flat_mesh = np.array([mesh.flatten() for mesh in meshes]).T
        default_params = " %d %d %d %d %s" % (
            value['defaults']['threshold'], value['defaults']['num_particles'], value['defaults']['M_t'],
            value['defaults']['M_f'], identifier)
        cmd_fmt = " ".join(["%.4f"] * flat_mesh.shape[1])
        with open(identifier + ".in", 'w') as f:
            for i in range(resolution ** 3):
                cmd_string = cmd_fmt % tuple(flat_mesh[i, :]) + default_params
                f.write(cmd_string + "\n")


if __name__ == '__main__':
    filename = 'params_single_channel.json'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    with open(filename, 'r') as param_file:
        parameter_sets = json.load(param_file)
    single_channel_exploration_set(parameter_sets)
    exit(0)
    for size in parameter_sets.keys():
        default_values = parameter_sets[size]['defaults']
        for i, param in enumerate(parameter_sets[size]['relations']):
            single_channel_param_set(param, default_values, "single_channel_data/param_file_%s_%d" % (size, i))
