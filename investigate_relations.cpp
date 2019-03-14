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

void test_constant_in_density() {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
     int num_steps = 20;
     int repeats = 25;
     double polarisation_ratio = 0.90;
     double final_time = 1E5;

     double gate_radius = 0.3;
     int gate_capacity = 2;
     double relation_constant = 1./200;
     for (int step=1;step < num_steps;step++) {
         int num_particles = 200 + step*10;
         double radius = std::sqrt(num_particles*relation_constant);
         printf("Number of particles: %d\nUrn Radius: %.2f\n",num_particles,radius);
         double pol_time = 0;
         for (int rep=0;rep<repeats;rep++) {
             Simulation sim = Simulation(num_particles,gate_radius);
             sim.left_gate_capacity = gate_capacity;
             sim.right_gate_capacity = gate_capacity;
             sim.circle_radius = radius;
             sim.setup();
             sim.start_evenly();
             int diff = 0;
             while (diff < num_particles*polarisation_ratio and sim.time < final_time) {
                 sim.update(0.0);
                 unsigned long it = sim.total_left.size() - 1;
                 diff = std::abs((int)sim.total_left.at(it) - (int)sim.total_right.at(it));
             }
             pol_time += sim.time;
         }
         pol_time /= repeats;
         printf("polarisation time: %.2f\n\n",pol_time);
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

double test_parameters(double gate_radius, int gate_capacity) {
    int repeats = 1000;
    double total_time = 0;
    for (int i = 0; i < repeats; i++) {
        total_time += get_thermalisation_time(gate_radius, gate_capacity);
    }
    return total_time / repeats;
}

int main(int argc, char *argv[]) {
    test_constant_in_density();
    return 0;
}

