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
    while (simulation.time < 100) {
        //simulation.print_status();
        simulation.update();
    }
    simulation.finish();
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
            run_domain();
            break;
        }
    }
    return 0;
}

