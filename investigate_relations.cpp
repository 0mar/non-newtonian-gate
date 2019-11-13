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
    results_file.close();
}

double get_chi(const unsigned long M_t, const unsigned long M_f, const double channel_length,
               const double channel_width, const double urn_radius, const int threshold, const int num_particles,
               const std::string &id) {
    const bool write_all_chi = true;
    double chi = 0;
    Simulation sim = Simulation(num_particles, channel_width, urn_radius, channel_length, threshold, threshold);
    sim.gate_is_flat = true;
    sim.distance_as_channel_length = true;
    sim.setup();
    std::random_device rd;
    std::mt19937 re(rd());
    std::uniform_real_distribution<double> unif(0.5, 1);
    std::ostringstream s;
    const double left_ratio = unif(re);
    try {
        sim.start(left_ratio);
    } catch (const std::invalid_argument &ex) {
        printf("Not running for bridge height %.2f and radius %.2f, returning 0\n", channel_width, urn_radius);
        return 0;
    }

    while (sim.measuring_times.size() < M_t) {
        sim.update(0.0);
    }
    if (write_all_chi) {
        s << sim.measuring_times.size() << "," << sim.get_mass_spread() << "," << channel_width << "," << urn_radius
          << "," << channel_length << "," << threshold;
    }
    const double weight = 1. / (double) (M_f - M_t);
    while (sim.measuring_times.size() < M_f) {
        sim.update(0.0);
        chi += weight * sim.get_mass_spread();
    }
    if (write_all_chi) {
        s << sim.measuring_times.size() << "," << sim.get_mass_spread() << std::endl;
        std::ofstream result_file(id + ".chi", std::ios::app);
        result_file << s.str();
        result_file.close();
    }
    sim.finish();
    return chi;
}


void matteo_relation_finder(int argc, char *argv[]) {
    const int num_arguments = 8;
    const int num_runs = 10;
    if (argc != num_arguments + 1) {
        std::cout << "Printing arguments: " << argc << std::endl;
        for (unsigned int i = 0; i < argc; i++) {
            std::cout << argv[i] << " ";
        }
        std::cout << std::endl;
        throw std::invalid_argument(
                "Please provide (in order) (1) channel length, (2) width, (3) urn radius, (4) threshold, "
                "(5) number of particles, (6) start point, (7) end point, (8) ID");
    }
    const double channel_length = std::stod(argv[1]);
    const double channel_width = std::stod(argv[2]);
    const double urn_radius = std::stod(argv[3]);
    const int threshold = std::stoi(argv[4]);
    const int num_particles = std::stoi(argv[5]);
    const int M_t = std::stoi(argv[6]);
    const int M_f = std::stoi(argv[7]);
    const std::string id = argv[8];
    double av_chi = 0;
    for (unsigned int i = 0; i < num_runs; i++) {
        av_chi += get_chi(M_t, M_f, channel_length, channel_width, urn_radius, threshold, num_particles, id) / num_runs;
    }
    std::ostringstream s;
    s << channel_length << "," << channel_width << "," << urn_radius << "," << threshold << "," << av_chi << std::endl;
    std::ofstream result_file(id + ".out", std::ios::app);
    result_file << s.str();
    result_file.close();

}

int main(int argc, char *argv[]) {
    matteo_relation_finder(argc, argv);
    return 0;
}

