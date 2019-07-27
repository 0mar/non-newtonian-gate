//
// Created by Omar Richardson on 24/10/2018.
//

#ifndef BLOTZ_SIMULATION_H // todo: fix guard name
#define BLOTZ_SIMULATION_H

#include <iostream>
#include <fstream>
#include <random>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>

struct Impact {
    unsigned long particle;
    double time;
    double stamp;

    bool operator<(const Impact &i) const {
        return time > i.time;
    }

    Impact(const unsigned long particle, const double time, const double stamp) : particle(particle), time(time),
                                                                                  stamp(stamp) {

    }
};

class Simulation {
public:
    Simulation(const int &num_particles, const double &gate_radius);

    // Important parameters
    const int num_particles;
    const double gate_radius;
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
    const unsigned long LEFT = 0;
    const unsigned long RIGHT = 1;
    bool explosion_direction_is_random;

    /**
     * Compute necessary parameters for the simulation, initialize data structures.
     * Run only once per simulation. Different runs require new setups and (therefore) new objects.
     */
    void setup();

    /**
     * Start the simulation. Initialize particles and times, and compute the next (first) impact.
     */
    void start();

    /**
     * Different start method (could possibly be merged).
     * As opposed to `start()`, this method puts as much particles left as it does right
     */
    void start_evenly();

    /**
     * Check if the point (x,y) is in the domain (collection for saying it is either in the left, the right or the gate.
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point
     * @return true if inside the domain, false otherwise.
     */
    bool is_in_domain(const double &x, const double &y); // todo: make this and other const

    /**
     * Check if the point (x,y) is in a circle on side `side`.
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point
     * @param side `LEFT` or `RIGHT`
     * @return `true` if point in circle, `false` otherwise
     */
    bool is_in_circle(const double &x, const double &y, const unsigned long &side);

    /**
     * Check if point is in the bridge
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point
     * @return `true` if point in bridge, `false` otherwise
     */
    bool is_in_bridge(const double &x, const double &y);

    /**
     * Compute the time it takes for a particle to reach the bridge.
     * If collision with bridge will not happen, return a time exceeding maximum collision time.
     * @param particle Particle index
     * @param normal_angle normal angle of the collision surface (output variable).
     * Only well-defined for not-exceeding-maximum-time collision time
     * @return time to next collision with bridge
     */
    double time_to_hit_bridge(const unsigned long &particle, double &normal_angle);

    /**
     * Computes the time it takes for a particle to reach the boundary of a reservoir
     * @param particle Particle index
     * @param center_x Center of the (circular) reservoir
     * @param normal_angle normal angle of the collision surface (output variable).
     * Only well-defined for not-exceeding-maximum-time collision time
     * @return time to next collision with circle
     */
    double time_to_hit_circle(const unsigned long &particle, const double &center_x, double &normal_angle);

    /**
     * Computes the time it takes for a particle to hit the gate.
     * No normal angle required because this collision does not reflect.
     * @param particle Particle index
     * @return time to next collision with gate
     */
    double time_to_hit_gate(const unsigned long &particle);

    /**
     * Time for a particle to hit the central vertical axis
     * @param particle Particle index
     * @return time to passing middle
     */
    double time_to_hit_middle(const unsigned long &particle);

    /**
     * Computes the next impact of a particle by finding the minimum impact time of all options.
     * @param particle Particle index
     */
    void compute_next_impact(const unsigned long &particle);

    double time;
    double last_written_time;
    std::vector<Impact> next_impact_times;
    std::vector<double> last_written_times;
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

    /**
     * Compute the current position of a particle, based on the `time` variable
     * Interpolates on the assumption that the particle has no collision from its last update until this time
     * (which should be a valid assumption) and extrapolates linearly for an exact position.
     * @param particle Particle index
     * @param x output variable for x position
     * @param y output variable for y position
     */
    void get_current_position(const unsigned long &particle, double &x, double &y);

    /**
     * Check if particle can enter gate. If gate is below threshold, enters the particle in the gate
     * If the gate exceeds the threshold, explodes the gate
     * @param particle Particle
     * @param direction Which gate is being accessed, LEFT or RIGHT.
     */
    void check_gate_admission(const unsigned long &particle, const unsigned long &direction);

    /**
     * Explode gate for a particle: give particle in the gate a reverse velocity.
     * @param particle Particle index
     * @param direction Which gate it explodes in, LEFT or RIGHT
     */
    void explode_gate(const unsigned long &particle, const unsigned long &direction);

    /**
     * Check if a point is in the LEFT or RIGHT gate.
     * @param x x-coordinate of point
     * @param y y-coordinate of point
     * @param direction Direction, LEFT or RIGHT
     * @return True if position in gate, false otherwise
     */
    bool is_in_gate(const double &x, const double &y, const unsigned long &direction);

    /**
     * Remove a particle from the gate.
     * @param particle particle index
     * @param direction LEFT or RIGHT
     */
    void check_gate_departure(const unsigned long &particle, const unsigned long &direction);

    /**
     * Compute the next collision, collide and update all particle positions.
     * Write to file if needed.
     * @param write_dt If positive, interpolate and write the positions of particles every `write_dt` to file.
     * If zero, don't write at all.
     */
    void update(const double &write_dt);

    /**
     * Retrieves which particle will collide next and at what time.
     * @param particle particle that will have the next collision
     * @param next_impact time at which the next collision will take place
     */
    void get_next_impact(unsigned long &particle, double &next_impact);

    /**
     * Store a time stamp and the number of particles left and right.
     */
    void measure();

    /**
     * Print the current status of the simulation to stdout
     */
    void print_status();

    /**
     * Write the positions to file at a certain time, which *should* be between the current time and the next collision.
     * However, this is not enforced in this method.
     * @param time Time at which positions would be interpolated.
     */
    void write_positions_to_file(const double &time);

    /**
     * Write all timestamps and number of positions to file.
     */
    void write_totals_to_file();

    /**
     * Compute the reflection angle based on an ingoing angle and the normal angle of the surface.
     * @param angle_in incoming angle of the particle
     * @param normal_angle normal angle of the surface
     * @return outgoing angle of the particle
     */
    double get_reflection_angle(const double &angle_in, const double &normal_angle);

    /**
     * Compute the angle a particle should have after a gate explosion (reverse/random)
     * @param particle Particle index
     * @return angle of the particle after exploding the gate
     */
    double get_retraction_angle(const unsigned long &particle);

    /**
     * Finish up simulation (write results, optional post-processing)
     */
    void finish();

private:
    /**
     * Compute the distance between the two reservoirs on bridge boundary height.
     * This is slightly larger than the distance between reservoir on bridge middle hight.
     */
    void couple_bridge();


    std::shared_ptr<std::random_device> rd;
    std::shared_ptr<std::mt19937> rng;
    std::shared_ptr<std::uniform_real_distribution<double>> unif_real;

};


#endif //BLOTZ_SIMULATION_H
