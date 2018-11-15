#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>
#include <fstream>

void single_particle_animation() {
    printf("Running the animation for a single particle\n");
    Simulation simulation = Simulation(1, 0);
    simulation.bridge_height = 0.5;
    simulation.setup();
    simulation.start();
    simulation.positions(0, 0) = -1.1;
    simulation.positions(0, 1) = 0.34;
    simulation.directions(0) = -2.3;
    simulation.compute_next_impact(0);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
}

void many_particle_animation() {
    printf("Running the animation for 200 particle\n");
    Simulation simulation = Simulation(200, 0);
    simulation.bridge_height = 0.5;
    simulation.setup();
    simulation.start();
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
}

void standard_simulation(double dt) {
    Simulation simulation = Simulation(100, 0.3);
    simulation.setup();
    simulation.start();
    simulation.write_positions_to_file(0);
    while (simulation.time < 500) {
        simulation.update(dt);
    }
    simulation.finish();
}

double get_thermalisation_time(double gate_radius, int gate_capacity) {
    Simulation simulation = Simulation(100, gate_radius);
    simulation.gate_capacity = gate_capacity;
    simulation.setup();
    simulation.start();
    simulation.write_positions_to_file(0);
    while (simulation.in_right.sum() < simulation.num_particles / 2) {
        simulation.update(0.05);
    }
    return simulation.time;
}

double test_parameters(double gate_radius, int gate_capacity) {
    int repeats = 10;
    double total_time = 0;
    for (int i = 0; i < repeats; i++) {
        total_time += get_thermalisation_time(gate_radius, gate_capacity);
    }
    return total_time / repeats;
}

int main(int argc, char *argv[]) {

    int mode = 0;
    if (argc == 2) {
        mode = std::stoi(argv[1]);
    } else if (argc == 3) {
        double gate_radius = std::stod(argv[1]);
        int gate_capacity = std::stoi(argv[1]);
        printf("Assuming a gate radius of %.2f and a gate capacity of %d\n", gate_radius, gate_capacity);
        test_parameters(gate_radius, gate_capacity);
    }
    switch (mode) {
        case 1: {
            single_particle_animation();
            break;
        }
        case 2: {
            many_particle_animation();
            break;
        }
        default: {
            std::cout << get_thermalisation_time(0.5, 6) << std::endl;
            break;
        }
    }
    return 0;
}

