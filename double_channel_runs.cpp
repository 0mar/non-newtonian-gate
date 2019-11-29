#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>
#include <fstream>
#include <ctime>
#include <sstream>


void get_chi(const double channel_width, const double channel_length, const int threshold, const double radius,
             const double second_width, const double second_length, const int num_particles, const double left_ratio,
             const unsigned long M_t, const unsigned long M_f, const std::string &id, double &av_chi, double &current) {
    Simulation sim = Simulation(num_particles, channel_width, radius, channel_length, threshold, threshold);
    sim.gate_is_flat = true;
    sim.distance_as_channel_length = true;
    sim.expected_collisions = M_f;
    sim.second_length = second_length;
    sim.second_width = second_width;
    sim.setup();
    std::ostringstream s;
    av_chi = 0;
    current = 0;
    try {
        sim.start(left_ratio);
    } catch (const std::invalid_argument &ex) {
        printf("Not running for bridge width %.2f and radius %.2f, returning 0\n", channel_width, radius);
    }

    while (sim.measuring_times.size() < M_t) {
        sim.update(0.0);
    }
    const double weight = 1. / (double) (M_f - M_t);
    const int count_offset = sim.first_channel_surplus;
    const double time_offset = sim.time;
    while (sim.measuring_times.size() < M_f) {
        sim.update(0.0);
        av_chi += weight * sim.get_mass_spread();
    }
    current = (sim.first_channel_surplus - count_offset) / (sim.time - time_offset);
}

//double get_chi_development(const unsigned long M_t, const unsigned long M_f, const double channel_length,
//                           const double channel_width, const double urn_radius, const int threshold,
//                           const double second_length, const double second_width, const int num_particles,
//                           const std::string &id) {
//
//    double chi = 0;
//    const int num_points = 500;
//    const unsigned long step_size = M_f / num_points;
//    Simulation sim = Simulation(num_particles, channel_width, urn_radius, channel_length, threshold, threshold);
//    sim.gate_is_flat = true;
//    sim.distance_as_channel_length = true;
//    sim.second_width = second_width; // todo: Choose width or width name
//    sim.second_length = second_length;
//    sim.setup();
//    std::random_device rd;
//    std::mt19937 re(rd());
//    std::uniform_real_distribution<double> unif(0.5, 1);
//    std::ostringstream s;
//    const double left_ratio = 0.25; // Not random atm
//    try {
//        sim.start(left_ratio);
//        sim.write_positions_to_file(0);
//    } catch (const std::invalid_argument &ex) {
//        printf("Not running for bridge width %.2f and radius %.2f, returning 0\n", channel_width, urn_radius);
//        return 0;
//    }
//    double dt = 0;
//    while (sim.measuring_times.size() < M_f) {
////        if (sim.measuring_times.size() ==M_t) {
////            dt = 0.025;
////            sim.last_written_time = sim.time;
////        }
////        if (sim.measuring_times.size() == M_t + 10000){
////            dt = 0;
////        }
//        sim.update(dt);
//        if (sim.measuring_times.size() % step_size == 0) {
//            s << sim.measuring_times.size() << "," << sim.time << "," << sim.first_channel_surplus << ","
//              << sim.second_channel_surplus << "," << sim.in_left << "," << std::fabs(sim.get_mass_spread())
//              << std::endl;
//        }
//    }
//    std::ofstream result_file(id + ".chi", std::ios::app);
//    result_file << s.str();
//    result_file.close();
//    return chi;
//}
// todo: Try a double loop

void mass_spread_and_current_for(int argc, char *argv[]) {
    const int num_arguments = 11;
    const int num_runs = 1;
    if (argc != num_arguments + 1) {
        std::cout << "Printing arguments: " << argc << std::endl;
        for (unsigned int i = 0; i < argc; i++) {
            std::cout << argv[i] << " ";
        }
        std::cout << std::endl;
        throw std::invalid_argument(
                "Please provide (in order) (1) channel width, (2) channel length, (3) threshold, (4) urn radius,"
                " (5) second channel width, (6) second channel length, (7) number of particles, (8) initial ratio,"
                " (9) transient time, (10) final time, (11) identifier");
    }
    const double channel_width = std::stod(argv[1]);
    const double channel_length = std::stod(argv[2]);
    const int threshold = std::stoi(argv[3]);
    const double radius = std::stod(argv[4]);
    const double second_width = std::stod(argv[5]);
    const double second_length = std::stod(argv[6]);
    const int num_particles = std::stoi(argv[7]);
    const double initial_ratio = std::stod(argv[8]);
    const int M_t = std::stoi(argv[9]);
    const int M_f = std::stoi(argv[10]);
    const std::string id = argv[11];
    double av_chi = 0;
    double current = 0;
    for (unsigned int i = 0; i < num_runs; i++) {
        get_chi(channel_length, channel_width, threshold, radius, second_length, second_width, num_particles,
                initial_ratio, M_t, M_f, id, av_chi, current);
    }
    std::ostringstream s;
    s << threshold << "," << second_width << "," << second_length << "," << initial_ratio << "," << av_chi
      << "," << current << std::endl;
    std::ofstream result_file(id + ".out", std::ios::app);
    result_file << s.str();
    result_file.close();
}

int main(int argc, char *argv[]) {
    mass_spread_and_current_for(argc, argv);
    return 0;
}

