#include <iostream>
#include "simulation.h"
#include <string>
#include <chrono>

void single_particle_animation() {
    printf("Running the animation for a single particle\n");
    Simulation simulation = Simulation(1, 0);
    simulation.bridge_height = 0.5;
    simulation.setup();
    simulation.start(1);
    simulation.x_pos.at(0) = -1.1;
    simulation.y_pos.at(0) = 0.34;
    simulation.directions.at(0) = -2.3;
    simulation.compute_next_impact(0);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
}

void many_particle_animation() {
    printf("Running the animation for 1000 particles, writing animation\n");
    Simulation simulation = Simulation(1000, 0.5);
    simulation.left_gate_capacity = 5;
    simulation.gate_is_flat = true;
    simulation.right_gate_capacity = 5;
    simulation.circle_distance = 0.5;
    simulation.circle_radius = 1;
    simulation.distance_as_channel_length = true;
    simulation.setup();
    simulation.start(0.25);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
    simulation.finish();
}

void time_test(const int num_times = 5) { // Starting point: 7 seconds, after keeping track of particle pos: 2.3
    printf("Running the animation for 10000 particles, timing\n");
    double average_time = 0;
    for (unsigned long i = 0; i < num_times; i++) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        Simulation simulation = Simulation(10000, 0.5);
        simulation.left_gate_capacity = 5;
        simulation.gate_is_flat = true;
        simulation.right_gate_capacity = 5;
        simulation.circle_distance = 0.5;
        simulation.circle_radius = 1;
        simulation.distance_as_channel_length = true;
        simulation.setup();
        simulation.start(0.25);
        while (simulation.measuring_times.size() < 1E6) {
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
        while (simulation.measuring_times.size() < 1E6) {
            if (simulation.measuring_times.size() % 100 == 0) {
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
            mass_spread_demo(20);
            break;
        }
        case 2: {
            time_test();
            break;
        }
        default: {
            many_particle_animation();
            break;
        }
    }
    return 0;
}

