//
// Created by Omar Richardson on 24/10/2018.
//

#ifndef BLOTZ_SIMULATION_H
#define BLOTZ_SIMULATION_H

#include <iostream>
#include <fstream>
#include <random>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>


class Simulation {
public:
    Simulation(int num_particles, double gate_radius);

    // Important parameters
    const int num_particles;
    double gate_radius;
    int left_gate_capacity;
    int right_gate_capacity;
    unsigned long in_left;
    unsigned long in_right;
    // Other parameters
    double circle_radius;
    double circle_distance;
    double bridge_height;
    // Computed quantities
    double left_center_x;
    double right_center_x;
    double max_path;
    double bridge_size;
    const int LEFT = 0;
    const int RIGHT = 1;

    void setup();

    void start();

    bool is_in_domain(double x, double y);

    bool is_in_left_circle(double x, double y);

    bool is_in_right_circle(double x, double y);

    bool is_in_bridge(double x, double y);

    double time_to_hit_bridge(unsigned long particle, double &normal_angle);

    double time_to_hit_circle(unsigned long particle, double center_x, double &normal_angle);

    double time_to_hit_gate(unsigned long particle);

    double time_to_hit_middle(unsigned long particle);

    double time;
    double last_written_time;
    std::vector<double> next_impact_times;
    std::vector<double> impact_times;
    std::vector<double> next_x_pos;
    std::vector<double> next_y_pos;
    std::vector<double> x_pos;
    std::vector<double> y_pos;
    std::vector<double> directions;
    std::vector<double> next_directions;
    std::vector<unsigned long> in_left_gate;
    std::vector<unsigned long> in_right_gate;

    std::vector<double> measuring_times;
    std::vector<unsigned long> total_left;
    std::vector<unsigned long> currently_in_left_gate;
    std::vector<unsigned long> currently_in_right_gate;
    std::vector<unsigned long> total_right;
    std::vector<std::vector<unsigned long>> gate_contents;
    std::vector<std::vector<unsigned long>> gate_arrays;
    std::vector<double> gate_capacities;

    void compute_next_impact(unsigned long particle);

    void get_current_position(unsigned long particle, double &x, double &y);

    void check_gate_admission(unsigned long particle, unsigned long direction);

    void explode_gate(unsigned long particle, unsigned long direction);

    bool is_in_gate(double x, double y, unsigned long direction);

    void check_gate_departure(unsigned long particle, unsigned long direction);

    void update(double write_dt);

    void measure();

    void print_status();

    void write_positions_to_file(double time);

    void write_totals_to_file();

    double get_reflection_angle(double angle_in, double normal_angle);

    double get_retraction_angle(unsigned long particle);

    void finish();

private:
    void couple_bridge();


    std::shared_ptr<std::random_device> rd;
    std::shared_ptr<std::mt19937> rng;
    std::shared_ptr<std::uniform_real_distribution<double>> unif_real;

};


#endif //BLOTZ_SIMULATION_H
