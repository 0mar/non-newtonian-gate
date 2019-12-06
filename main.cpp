#include <iostream>
#include "simulation.h"
#include <string>
#include <chrono>

/**
 * This file contains recipes for the most common executables for the two-chamber system.
 * The recipes themselves are pretty self-explanatory.
 */

void polarisation_demo() {
    printf("Creating the animation for a polarising system\n");
    Simulation simulation = Simulation(300, 0.2);
    simulation.left_gate_capacity = 2;
    simulation.gate_is_flat = true;
    simulation.right_gate_capacity = 2;
    simulation.circle_distance = 0.5;
    simulation.circle_radius = 0.5;
    simulation.distance_as_channel_length = true;
    simulation.setup();
    simulation.start(0.5);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (std::fabs(simulation.get_mass_spread()) < 0.9 and simulation.time < 1000) {
        simulation.update(dt);
    }
    printf("Polarised in %.2f seconds, mass spread of %.2f\n", simulation.time, simulation.get_mass_spread());
}

void double_channel_demo() {
    printf("Creating the animation for 1000 particles\n");
    Simulation simulation = Simulation(1000, 0.5);
    simulation.left_gate_capacity = 7;
    simulation.gate_is_flat = true;
    simulation.right_gate_capacity = 7;
    simulation.circle_distance = 0.5;
    simulation.circle_radius = 1;
    simulation.second_length = 1;
    simulation.second_width = 0.3;
    simulation.distance_as_channel_length = true;
    simulation.setup();
    simulation.start(0.75);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
}

void many_particle_animation() {
    printf("Creating the animation for 1000 particles\n");
    Simulation simulation = Simulation(100, 0.5);
    simulation.left_gate_capacity = 3;
    simulation.gate_is_flat = true;
    simulation.right_gate_capacity = 3;
    simulation.circle_distance = 0.5;
    simulation.circle_radius = 1;
    simulation.second_length = 1;
    simulation.second_width = 0.3;
    simulation.distance_as_channel_length = true;
    simulation.setup();
    simulation.start(0);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
}

void time_test(const int num_times = 5) { // Starting point: 7 seconds, after keeping track of particle pos: 2.3
    printf("Running the animation for 10000 particles, timing\n");
    double average_time = 0;
    for (unsigned long i = 0; i < num_times; i++) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        Simulation simulation = Simulation(10000, 0.5);
        simulation.left_gate_capacity = 3;
        simulation.gate_is_flat = true;
        simulation.right_gate_capacity = 3;
        simulation.circle_distance = 0.5;
        simulation.circle_radius = 1;
        simulation.distance_as_channel_length = true;
        simulation.setup();
        simulation.start(0.25);
        while (simulation.num_collisions < (unsigned long) 1E6) {
            simulation.update(0);
        }
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_span = t2 - t1;
        average_time += time_span.count() / num_times;
    }
    printf("Time took: %.3f ms\n", average_time);
}

void mass_spread_demo(const int &num_res) {
    printf("Running the simulation for 1000 particles (million collisions), writing mass spread\n");
    std::ofstream results_file;
    std::ostringstream s;
    for (unsigned long j = 0; j < num_res; j++) {
        Simulation simulation = Simulation(1000, 0.5);
        simulation.left_gate_capacity = 4;
        simulation.gate_is_flat = true;
        simulation.right_gate_capacity = 4;
        simulation.circle_distance = 0.5;
        simulation.circle_radius = 0.7;
        simulation.distance_as_channel_length = true;
        simulation.write_bounces = false;
        simulation.setup();
        simulation.start(0.7);
        simulation.write_positions_to_file(0);
        while (simulation.num_collisions < (unsigned long) 1E6) {
            if (simulation.num_collisions % 100 == 0) {
                s << simulation.get_mass_spread() << ",";
            }
            simulation.update(0);
        }
        s << simulation.get_mass_spread();
        s << std::endl;
    }
    results_file.open("mass_spread_over_time.txt");
    results_file << s.str();
    results_file.close();
}

int main(int argc, char *argv[]) {
    int mode = 0;
    if (argc == 2) {
        mode = std::stoi(argv[1]);
    }
    switch (mode) {
        case 1: {
            double_channel_demo();
            break;
        }
        case 2: {
            many_particle_animation();
            break;
        }
        default: {
            polarisation_demo();
            break;
        }
    }
    return 0;
}

