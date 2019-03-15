#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>
#include <fstream>

template<class T>
void write_results(std::string &id, std::vector<T> &data) {
    std::ofstream results_file;
    results_file.open(id + ".txt");
    for (T datum: data) {
        results_file << datum << "\t";
    }
    results_file << std::endl;
}

int get_critical_number_of_particles(double radius, int capacity, double gate_height, double gate_length,
                                     int lower_bound = 1000, int upper_bound = 7000, int guess = 3500) {
    double final_time = 1E5;
    double polarisation_ratio = 0.95;
    double gate_radius = gate_length / 2;
    int repeats = 3;
    bool found_lower_bound = false;
    bool found_upper_bound = false;
    int num_particles = guess;
    while (not (found_lower_bound and found_upper_bound)) {
        printf("Testing %d particles\n", num_particles);
        int num_of_polarizations=0;
        for (int rep=0;rep < repeats;rep++) {
            Simulation sim = Simulation(num_particles,gate_radius);
            sim.left_gate_capacity = capacity;
            sim.right_gate_capacity = capacity;
            sim.circle_radius = radius;
            sim.bridge_height = gate_height;
            sim.circle_distance = gate_length; // todo: Remember: approximative.
            sim.setup();
            sim.start_evenly();
            int diff = 0;
            while (diff < num_particles*polarisation_ratio and sim.time < final_time) {
                sim.update(0.0);
                unsigned long it = sim.total_left.size() - 1;
                diff = std::abs((int)sim.total_left.at(it) - (int)sim.total_right.at(it));
            }
            bool has_polarized = (sim.time < final_time);
            if (has_polarized) {
                num_of_polarizations++;
                printf("has polarized in %.2f\n", sim.time);
            } else {
                printf("not polarized in %.2f\n", final_time);
            }
        }
            if (num_of_polarizations==0) {
                lower_bound = num_particles;
                num_particles = (lower_bound + upper_bound) / 2;
                if (num_particles - num_of_polarizations < 10) {
                    printf("Very close to lower bound and still not polarizing; reducing lower bound\n");
                    lower_bound -= 50;
                }
            } else if (num_of_polarizations==repeats){
                upper_bound = num_particles;
                num_particles = (lower_bound + upper_bound) / 2;
                if (upper_bound - num_particles < 10) {
                    printf("Very close to upper bound and still not polarizing: reducing upper bound\n");
                    upper_bound += 50;
                }
            } else {
                printf("Interval: (%d,%d), crit: %d\n", lower_bound, upper_bound, num_particles);
                found_lower_bound = true;
                found_upper_bound = true;
            }
        }
    return num_particles;
}

void test_constant_in_density(int capacity) {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
     int num_steps = 20;
    double height = 0.3;
    double length = 0.6;
    printf("Testing constant in density: Capacity:%d/\n", capacity);
    int lb = 0;
    int ub = 2000; // Link these three to your initial guess by some estimate; count on continuity
    int guess = 1000;
     for (int step=0;step < num_steps;step++) {
         double radius = 2 + step*0.2;
         int crit = get_critical_number_of_particles(radius, capacity, height, length, lb, ub, guess);
         printf("o: Radius:%.2f/Critical Number:%d/\n", radius, crit);
         guess = (int) (crit * 1.5);
         lb = (int) (crit * 0.9);
         ub = (int) (crit * 1.6);
     }
}

void test_linear_in_capacity(double density) {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
    int num_steps = 20;
    double radius = 2.;
    double height = 0.3;
    double length = 0.6;
    printf("Testing linear in capacity\n");
    int lb = 0;
    int ub = 2000; // Link these three to your initial guess by some estimate; count on continuity
    int guess = 1000;
    for (int step = 0; step < num_steps; step++) {
        int capacity = 2 + step * 2;
        int crit = get_critical_number_of_particles(radius, capacity, height, length, lb, ub, guess);
        printf("o: Capacity:%d/Critical Number:%d/\n", capacity, crit);
        guess = (int) (crit * 1.3);
        lb = (int) (crit * 0.9);
        ub = (int) (crit * 1.6);
    }
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
int main(int argc, char *argv[]) {
    int mode = 0;
    if (argc == 2) {
        mode = std::stoi(argv[1]);
    }
    switch (mode) {
        case 1: {
            test_constant_in_density(2);
            break;
        }
        case 2: {
            test_constant_in_density(20);
            break;
        }
        case 3: {
            test_linear_in_capacity(0.5);
            break;
        }
        default: {
            printf("Please choose option (0-2)\n");
            break;
        }
    }
    return 0;
}

