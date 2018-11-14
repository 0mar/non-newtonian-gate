#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>
#include <fstream>

void write_results(std::string &id, std::vector<double> &x, std::vector<double> &y) {
    if (x.size() != y.size()) std::cout << "Writing fails..." << std::endl;
    std::ofstream results_file;
    results_file.open(id + ".txt");
    results_file << "[";
    for (unsigned long i = 0; i < x.size(); i++) {
        results_file << "[" << x.at(i) << ", " << y.at(i) << "]," << std::endl;
    }
    results_file << "]";
}

void run_domain() {
    Simulation simulation = Simulation(1000, 0.3);
    simulation.setup();
    simulation.start();
    simulation.write_to_file(false);
    while (simulation.time < 50) {
        //simulation.print_status();
        simulation.update(0);
        simulation.write_to_file(false);
    }
    simulation.finish();
}

void write_domain(double dt) {
    Simulation simulation = Simulation(100, 0.3);
    simulation.setup();
    simulation.start();
    simulation.write_to_file(true);
    int iteration = 0;
    while (simulation.time < 500) {
        iteration++;
        while (simulation.time < iteration * dt) {
            simulation.update(iteration * dt);
        }
    }
    //simulation.finish();
}


int main(int argc, char *argv[]) {
    int mode = 0;
    if (argc == 2) {
        mode = std::stoi(argv[1]);
    } else {
    }
    switch (mode) {
        case 1: {
            run_domain();
            break;
        }
        default: {
            write_domain(0.03);
            break;
        }
    }
    return 0;
}

