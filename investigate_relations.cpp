#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>
#include <fstream>
#include <ctime>
#include <sstream>

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
                                     int lower_bound = 1000, int upper_bound = 7000, int guess = 3500,
                                     bool high_res = true) {
    /*
     * Find the critical number of particles for any set of parameters.
     * This is a binary search method that from some (guess) lower and upper bound tries to derive an estimate
     * for the critical number of particles: where the situation sometimes thermalizes and sometimes does not.
     *
     * Because this is a chaotic system and all depends on the initial positions,
     * not observing critical behaviour does not mean that it is not there. To mediate for this, we try every guess three times.
     *
     * another problem is that the simulations take time and not the entire search space can feasibly be explored.
     * To resolve this, we take an indicative upper and lower bound and move them up when we cannot find the behaviour we are looking for.
     * This is done by first assuming no real upper bound until we observe antithermalisation, and assuming the real lower bound is zero.
     */
    double final_time = 10E4;
    double polarisation_ratio = 0.95;
    double gate_radius = gate_length / 2;
    // How often will we repeat the experiment?
    int repeats = 3;
    // How many deviant outcomes will we tolerate?
    int num_deviant = 0;
    if (high_res) {
        repeats = 7;
        num_deviant = 1;
    }
    bool found = false;
    bool found_upper_bound = false;
    int num_particles = guess;
    while (not found) {
        printf("Testing %d particles\n", num_particles);
        int num_of_polarizations = 0;
        for (int rep = 0; rep < repeats; rep++) {
            Simulation sim = Simulation(num_particles, gate_radius);
            sim.left_gate_capacity = capacity;
            sim.right_gate_capacity = capacity;
            sim.circle_radius = radius;
            sim.bridge_height = gate_height;
            sim.circle_distance = gate_length; // todo: Remember: approximative.
            sim.setup();
            sim.start_evenly();
            int diff = 0;
            while (diff < num_particles * polarisation_ratio and sim.time < final_time) {
                sim.update(0.0);
                unsigned long it = sim.total_left.size() - 1;
                diff = std::abs((int) sim.total_left.at(it) - (int) sim.total_right.at(it));
            }
            bool has_polarized = (sim.time < final_time);
            if (has_polarized) {
                num_of_polarizations++;
                printf("has polarized in %.2f\n", sim.time);
            } else {
                printf("not polarized in %.2f\n", final_time);
            }
        }
        if (num_of_polarizations <= num_deviant) {
            if (not found_upper_bound) {
                printf("Doubling upper bound\n");
                lower_bound = num_particles;
                num_particles = upper_bound;
                upper_bound = num_particles * 2;
            } else {
                lower_bound = num_particles;
                num_particles = (lower_bound + upper_bound) / 2;
                if (upper_bound - num_particles < 5) {
                    printf("Very close to upper bound and still not polarizing: increasing upper bound\n");
                    upper_bound += 10;
                    num_particles = (lower_bound + upper_bound) / 2;
                }
            }
        } else if (num_of_polarizations >= repeats - num_deviant) {
            upper_bound = num_particles;
            num_particles = (lower_bound + upper_bound) / 2;
            found_upper_bound = true;
            if (num_particles - lower_bound < 5) {
                printf("Very close to lower bound and still not polarizing; reducing lower bound\n");
                if (lower_bound < 20) {
                    printf("Lower bound of %d not low enough. No meaningful convergence, exiting\n", lower_bound);
                    std::exit(200);
                } else {
                    lower_bound -= 20;
                }
                num_particles = (lower_bound + upper_bound) / 2;
            }
        } else {
            printf("Interval: (%d,%d), crit: %d\n", lower_bound, upper_bound, num_particles);
            found = true;
        }
    }
    return num_particles;
}

void test_constant_in_density(int capacity, int guess) {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
    std::ofstream result_file("result_run" + std::to_string(std::time(nullptr)) + ".txt");
    int num_steps = 10;
    double height = 0.3;
    double length = 0.6;
    result_file << "number of steps: " << num_steps << "\theight: " << height << "\tlength " << length << "\tcapacity: "
                << capacity << std::endl;
    printf("Testing constant in density: Capacity:%d/\n", capacity);
    int lb = 0;
    int ub = guess * 2; // Link these three to your initial guess by some estimate; count on continuity
    for (int step = 0; step < num_steps; step++) {
        double radius = .35 + step * 0.05;
        int crit = get_critical_number_of_particles(radius, capacity, height, length, lb, ub, guess);
        printf("o: Radius:%.2f/Critical Number:%d/\n", radius, crit);
        result_file << "radius " << radius << " crit " << crit << std::endl;
        guess = (int) (crit * 1.5);
        lb = (int) (crit * 0.9);
        ub = (int) (crit * 1.6);
    }
    result_file.close();
}

void test_linear_in_capacity(double density) {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
    std::ofstream result_file("result_run" + std::to_string(std::time(nullptr)) + ".txt");
    int num_steps = 20;
    double radius = 1;
    double height = 0.3;
    double length = 0.6;
    result_file << "number of steps: " << num_steps << "\theight: " << height << "\tlength " << length << "\tradius: "
                << radius << std::endl;
    printf("Testing linear in capacity\n");
    int lb = 0;
    int ub = 600; // Link these three to your initial guess by some estimate; count on continuity
    int guess = 400;
    for (int step = 0; step < num_steps; step++) {
        int capacity = 2 + step;
        radius -= 0.05;
        int crit = get_critical_number_of_particles(radius, capacity, height, length, lb, ub, guess);
        printf("o: Capacity:%d/Critical Number:%d/\n", capacity, crit);
        result_file << "capacity " << capacity << " crit " << crit << std::endl;
        guess = (int) (crit * 1.3);
        lb = (int) (crit * 0.9);
        ub = (int) (crit * 1.6);
    }
    result_file.close();
}

void test_inverse_in_height() {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
    std::ofstream result_file("result_run" + std::to_string(std::time(nullptr)) + ".txt");
    int num_steps = 40;
    double radius = 1;
    double length = 0.6;
    int capacity = 2;
    result_file << "number of steps: " << num_steps << "\tcapacity: " << capacity << "\tlength " << length
                << "\tradius: "
                << radius << std::endl;
    printf("Testing inverse in height\n");
    int lb = 0;
    int ub = 1000; // Link these three to your initial guess by some estimate; count on continuity
    int guess = 500;
    for (int step = 0; step < num_steps; step++) {
        double height = 0.2 + step * 0.01;
        int crit = get_critical_number_of_particles(radius, capacity, height, length, lb, ub, guess);
        printf("o: Height:%.2f/Critical Number:%d/\n", height, crit);
        result_file << "Height " << height << " crit " << crit << std::endl;
        guess = (int) (crit * 1.3);
        lb = (int) (crit * 0.9);
        ub = (int) (crit * 1.6);
    }
    result_file.close();
}

void test_inverse_in_length() {
    /**
     * This method tests if we get constant thermalisation time for number of particles scaling with the density.
     * This method must be tested in the critical case because otherwise there is no thermalisation happening.
     */
    std::ofstream result_file("result_run" + std::to_string(std::time(nullptr)) + ".txt");
    int num_steps = 20;
    double radius = 1;
    double height = 0.3;
    int capacity = 2;
    result_file << "number of steps: " << num_steps << "\tcapacity: " << capacity << "\theight " << height
                << "\tradius: "
                << radius << std::endl;
    printf("Testing inverse in length\n");
    int lb = 0;
    int ub = 2000; // Link these three to your initial guess by some estimate; count on continuity
    int guess = 1000;
    for (int step = 0; step < num_steps; step++) {
        double length = 0.3 + step * 0.1;
        int crit = get_critical_number_of_particles(radius, capacity, height, length, lb, ub, guess);
        printf("o: Length:%.2f/Critical Number:%d/\n", length, crit);
        result_file << "length " << length << " crit " << crit << std::endl;
        guess = (int) (crit * 1.3);
        lb = (int) (crit * 0.9);
        ub = (int) (crit * 1.6);
    }
    result_file.close();
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

double
get_chi(const unsigned long M_t, const unsigned long M_f, const double channel_length, const double channel_width,
        const double urn_radius, const int threshold) {
    double chi = 0;
    const int num_particles = 1000;
    Simulation sim = Simulation(num_particles, 1, 1, channel_length, channel_width / 2, threshold, threshold);
    sim.setup();
    sim.start();
    while (sim.measuring_times.size() < M_t) {
        sim.update(0.0);
    }
    const double weight = 1. / (double) (M_f - M_t);
    while (sim.measuring_times.size() < M_f) {
        sim.update(0.0);
        chi += weight * std::fabs(sim.total_left.back() - sim.total_right.back()) / sim.num_particles;
    }
    return chi;
}

void omar_relation_finder(int argc, char *argv[]) {
    int mode = 0;
    if (argc == 2) {
        mode = std::stoi(argv[1]);
    }
    switch (mode) {
        case 1: {
            test_constant_in_density(2, 1000);
            break;
        }
        case 2: {
            test_constant_in_density(20, 750);
            break;
        }
        case 3: {
            test_linear_in_capacity(0.5);
            break;
        }
        case 4: {
            test_inverse_in_height();
            break;
        }
        case 5: {
            test_inverse_in_length();
            break;
        }
        default: {
            printf("Please choose option (1-5)\n");
            break;
        }
    }
}

void matteo_relation_finder(int argc, char *argv[]) {
    const int num_arguments = 5;
    const int num_runs = 4;
    if (argc != num_arguments + 1) {
        std::cout << "Printing arguments: " << argc << std::endl;
        for (unsigned int i = 0; i < argc; i++) {
            std::cout << argv[i] << " ";
        }
        std::cout << std::endl;
        throw std::invalid_argument("Please provide (in order) channel length, width, urn radius, threshold and ID");
    }
    const double channel_length = std::stod(argv[1]);
    const double channel_width = std::stod(argv[2]);
    const double urn_radius = std::stod(argv[3]);
    const int threshold = std::stoi(argv[4]);
    const std::string id = argv[5];
    double tot_chi = 0;
    const unsigned long M_t = 1E5; // Number of hits until we start measuring
    const unsigned long M_f = 1.5E5; // Number of hits until we stop measuring;
    for (unsigned int i = 0; i < num_runs; i++) {
        tot_chi += get_chi(M_t, M_f, channel_length, channel_width, urn_radius, threshold) / num_runs;
    }
    std::ostringstream s;
    s << channel_length << "," << channel_width << "," << urn_radius << "," << threshold << "," << tot_chi << std::endl;
    std::ofstream result_file(id + ".out", std::ios::app);
    result_file << s.str();
    result_file.close();

}

int main(int argc, char *argv[]) {
    matteo_relation_finder(argc, argv);
    return 0;
}

