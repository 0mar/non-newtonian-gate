#include <iostream>
#include <memory>
#include "simulation.h"
#include <string>

/**
 * This file contains an executable that can be used for efficient parameter regime explorations
 * of the two-chamber dynamics
 * It takes command line parameters that define the simulation, allowing for simple batch running
 * which allows for efficient parallel execution
 */

/**
 * Obtain the mass spread as a function of the parameters below, averaged from the transient time to the final time.
 * For a definition of the mass spread, see simulation.h.
 *
 * @param M_t Transient time, measured in number of collisions
 * @param M_f Final time, measured in number of collisions
 * @param channel_length Length of the channel
 * @param channel_width Width of the channel
 * @param urn_radius Radius of the chamber
 * @param threshold Number of particles that can at the same time in the channel
 * @param num_particles Number of particles in the system
 * @param id File identifier to write auxiliary results to.
 * @return Average mass spread
 */

double get_chi(const unsigned long M_t, const unsigned long M_f, const double channel_length,
               const double channel_width, const double urn_radius, const int threshold, const int num_particles,
               const std::string &id, const double &init_chi) {
    double chi = 0;
    Simulation sim = Simulation(num_particles, channel_width, urn_radius, channel_length, threshold, threshold);
    sim.gate_is_flat = true;
    sim.distance_as_channel_length = true;
    sim.expected_collisions = M_f;
    const double left_ratio = (init_chi+1)/2;

    std::ostringstream s;
    try {
        sim.setup();
    } catch (const std::invalid_argument &ex) {
        printf("Not running for bridge width %.2f and radius %.2f, returning 0\n", channel_width, urn_radius);
        return 0;
    }
    sim.start(left_ratio);
    while (sim.num_collisions < M_t) {
        sim.update(0.0);
    }
    const double weight = 1. / (double) (M_f - M_t);
    while (sim.num_collisions < M_f) {
        sim.update(0.0);
        chi += weight * sim.get_mass_spread();
    }
    chi = std::fabs(chi);
    return chi;
}

/**
 * Obtain the mass spread as a function of the parameters below, and write its value 500 times during the proces .
 * For a definition of the mass spread, see simulation.h.
 *
 * @param M_t Transient time, measured in number of collisions
 * @param M_f Final time, measured in number of collisions
 * @param channel_length Length of the channel
 * @param channel_width Width of the channel
 * @param urn_radius Radius of the chamber
 * @param threshold Number of particles that can at the same time in the channel
 * @param num_particles Number of particles in the system
 * @param id File identifier to write auxiliary results to.
 * @return Average mass spread
 */
double get_chi_development(const unsigned long M_t, const unsigned long M_f, const double channel_length,
                           const double channel_width, const double urn_radius, const int threshold,
                           const int num_particles, const std::string &id, const double &init_chi) {

    double chi = 0;
    const int num_points = 500;
    const unsigned long step_size = M_f / num_points;
    Simulation sim = Simulation(num_particles, channel_width, urn_radius, channel_length, threshold, threshold);
    sim.gate_is_flat = true;
    sim.distance_as_channel_length = true;
    sim.expected_collisions = M_f;
    sim.setup();
    std::random_device rd;
    std::mt19937 re(rd());
    std::uniform_real_distribution<double> unif(0.5, 1);
    std::ostringstream s;
    const double left_ratio = (init_chi+1)/2;
    try {
        sim.start(left_ratio);
        sim.write_positions_to_file(0);
    } catch (const std::invalid_argument &ex) {
        printf("Not running for bridge width %.2f and radius %.2f, returning 0\n", channel_width, urn_radius);
        return 0;
    }
    double dt = 0;
    while (sim.num_collisions < M_f) {
        sim.update(dt);
        if (sim.num_collisions % step_size == 0) {
            s << sim.num_collisions << "," << sim.time << "," << std::fabs(sim.get_mass_spread())
              << std::endl;
        }
    }
    std::ofstream result_file(id + ".chi", std::ios::app);
    result_file << s.str();
    result_file.close();
    return chi;
}

/**
 * Compute the time it takes for a fully polarized and fully equilibriated system to reach the same state.
 * This time can then be used as an indication for the thermilisation time of the system.
 *
 * @param channel_length
 * @param channel_width
 * @param urn_radius
 * @param threshold
 * @param num_particles
 * @param id
 * @param chi_diff Tolerance for difference in mass spreads
 * @return Time for the two systems to reach the same mass spread
 */
unsigned long find_sandwich_time(const double channel_length,
                                 const double channel_width, const double urn_radius, const int threshold,
                                 const int num_particles, const std::string &id, double &chi_diff) {
    double chi = 0;
    const int num_points = 500;
    const unsigned long M_t = 1E5;
    const unsigned long M_f = 1E8;
    const unsigned long step_size = M_f / num_points;
    Simulation sim_top = Simulation(num_particles, channel_width, urn_radius, channel_length, threshold, threshold);
    Simulation sim_bottom = Simulation(num_particles, channel_width, urn_radius, channel_length, threshold, threshold);
    sim_top.gate_is_flat = true;
    sim_bottom.gate_is_flat = true;
    sim_top.distance_as_channel_length = true;
    sim_bottom.distance_as_channel_length = true;
    sim_top.setup();
    sim_bottom.setup();
    std::ostringstream s;
    try {
        sim_top.start(0.5);
        sim_bottom.start(1);
    } catch (const std::invalid_argument &ex) {
        printf("Not running for bridge width %.2f and radius %.2f, returning 0\n", channel_width, urn_radius);
        return 0;
    }
    chi_diff = 1;
    const double eps = 1E-3;
    while (sim_top.num_collisions < M_t or (chi_diff > eps and sim_top.num_collisions < M_f)) {
        sim_top.update(0.0);
        sim_bottom.update(0.0);
        if (sim_top.num_collisions % step_size == 0) {
            s << sim_top.num_collisions << "," << std::fabs(sim_top.get_mass_spread()) << ", "
              << std::fabs(sim_bottom.get_mass_spread()) << std::endl;
        }
        chi_diff = std::fabs(std::fabs(sim_top.get_mass_spread()) - std::fabs(sim_bottom.get_mass_spread()));
    }
    std::ofstream result_file(id + ".chi", std::ios::app);
    result_file << s.str();
    result_file.close();
    return sim_top.num_collisions;
}

/**
 * Find the thermalisation times for a specific set of parameters. Executable, takes command line values.
 *
 * @param argc number of command line arguments
 * @param argv list of command line arguments
 */
void find_relation_time(int argc, char *argv[]) {
    const int num_arguments = 6;
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
    const std::string id = argv[6];
    double chi_diff = 0;
    unsigned long convergence_time = find_sandwich_time(channel_length, channel_width, urn_radius, threshold,
                                                        num_particles, id, chi_diff);
    std::ostringstream s;
    s << channel_length << "," << channel_width << "," << urn_radius << "," << threshold << "," << convergence_time
      << "," << chi_diff << std::endl;
    std::ofstream result_file(id + ".out", std::ios::app);
    result_file << s.str();
    result_file.close();
}

/**
 * Find the average mass spread in the thermalised state for a specific set of parameters.
 * Executable, takes command line values
 *
 * @param argc number of command line arguments. Must be 8
 * @param argv List of command line arguments. Must be in specific order.
 */
void average_mass_spread_for(int argc, char *argv[]) {
    const int num_arguments = 8;
    const int num_runs = 1;
    const int num_init_chis = 5;
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
    for (unsigned int j=0; j < num_init_chis; j++) {
        const double init_chi = j*0.25;
        for (unsigned int i = 0; i < num_runs; i++) {
            av_chi += get_chi(M_t, M_f, channel_length, channel_width, urn_radius, threshold, num_particles, id, init_chi)/num_runs;
        }
    std::ostringstream s;
    s << channel_length << "," << channel_width << "," << urn_radius << "," << threshold << "," << init_chi << "," << av_chi
      << std::endl;
    std::ofstream result_file(id + ".out", std::ios::app);
    result_file << s.str();
    result_file.close();
    }
}

int main(int argc, char *argv[]) {
    average_mass_spread_for(argc, argv);
    return 0;
}

