#include <iostream>
#include "simulation.h"
#include <string>

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
    printf("Running the animation for 200 particles, trying det exp\n");
    Simulation simulation = Simulation(1000, 0.5);
    simulation.left_gate_capacity = 4;
    simulation.gate_is_flat = true;
    simulation.right_gate_capacity = 4;
    simulation.circle_distance = 0.5;
    simulation.circle_radius = 0.7;
    simulation.distance_as_channel_length = true;

    simulation.setup();
    simulation.start(0.7);
    simulation.write_positions_to_file(0);
    double dt = 0.025;
    while (simulation.time < 100) {
        simulation.update(dt);
    }
    simulation.finish();
}

void standard_simulation(double dt) {
    Simulation simulation = Simulation(100, 0.3);
    simulation.setup();
    simulation.start(1);
    simulation.write_positions_to_file(0);
    while (simulation.time < 500) {
        simulation.update(dt);
    }
    simulation.finish();
}

void mass_spread_demo(const int &num_res) {
    printf("Running the animation for 1000 particles, trying det exp\n");

    std::ofstream results_file;
    results_file.open("chis.txt");

    std::vector<std::vector<double>> all_chis;
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
        double chi = 0;
        std::vector<double> chis;
        while (simulation.total_left.size() < 1E6) {
            simulation.update(0);
            chi = std::fabs(1. * simulation.total_left.back() - 1. * simulation.total_right.back()) /
                  simulation.num_particles;
            chis.push_back(chi);
        }
        all_chis.push_back(chis);
    }
    for (unsigned long el = 0; el < all_chis[0].size(); el++) {
        for (auto chis:all_chis) {
            results_file << chis.at(el) << ",";
        }
        results_file << std::endl;
    }
    results_file << std::endl;
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
            single_particle_animation();
            break;
        }
        default: {
            many_particle_animation();
            break;
        }
    }
    return 0;
}

