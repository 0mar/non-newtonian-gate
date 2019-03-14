#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>
#include <fstream>


void write_results(std::string &id, std::vector<double> &data) {
    std::ofstream results_file;
    results_file.open(id + ".txt");
    for (double datum: data) {
        results_file << datum << "\t";
    }
    results_file << std::endl;
}

void single_particle_animation() {
    printf("Running the animation for a single particle\n");
    Simulation simulation = Simulation(1, 0);
    simulation.bridge_height = 0.5;
    simulation.setup();
    simulation.start();
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
    printf("Running the animation for 200 particles, trying det exp\n");
    Simulation simulation = Simulation(50, 0.7);
    simulation.left_gate_capacity = 2;
    simulation.right_gate_capacity = 2;
    simulation.bridge_height = 0.5;
    simulation.explosion_direction_is_random = false;
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
    simulation.left_gate_capacity = gate_capacity;
    simulation.right_gate_capacity = gate_capacity;
    simulation.setup();
    simulation.start();
    // simulation.write_positions_to_file(0);
    while (simulation.total_right.at(simulation.total_right.size() - 1) < simulation.num_particles / 2 and
           simulation.time < 1E5) {
        simulation.update(0.0);
    }
    return simulation.time;
}

void converge_thermalisation(double gate_radius, int gate_capacity, unsigned long repeats) {
    auto results = std::vector<double>();
    results.reserve(repeats);
    for (unsigned long i = 0; i < repeats; i++) {
        results.push_back(get_thermalisation_time(gate_radius, gate_capacity));
        std::cout << i << std::endl;
    }
    std::string name = "therms";
    write_results(name, results);
}

double test_parameters(double gate_radius, int gate_capacity) {
    int repeats = 1000;
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
        int gate_capacity = std::stoi(argv[2]);
        //printf("Assuming a gate radius of %.2f and a gate capacity of %d\n", gate_radius, gate_capacity);
        std::cout << test_parameters(gate_radius, gate_capacity) << std::endl;
        return 0;
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
            converge_thermalisation(0.5, 5, 100000);
            break;
        }
    }
    return 0;
}

