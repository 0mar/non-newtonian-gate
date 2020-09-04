//
// Created by Omar Richardson on 24/10/2018.
//

#ifndef TERRIER_SIMULATION_H
#define TERRIER_SIMULATION_H

#include <iostream>
#include <fstream>
#include <random>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <ctime>
#include <sstream>
#include <unistd.h>
#include <numeric>

class Simulation {
public:
    /**
     * Create a new Simulation with given parameters.
     * Note that before `setup` has been called, the parameters can freely be changed,
     * allowing for some flexibility.
     * A lot of default parameters are present, the only ones of real relevance are `second_width`, which defaults to
     * 0 and then ignores the second channel and`flat_gate` which reinterprets the parameters to change the geometry.

     * @param num_particles Number of particles
     * @param bridge_width Width of the central channel
     * @param circle_radius Radius of the chamber (circle)
     * @param circle_distance Distance between the chambers
     * OR length of the bridge, based on `distance_as_channel_length`
     * @param left_gate_capacity Maximum number of particles allowed to enter the left channel
     * @param right_gate_capacity Maximum number of particles allowed to enter the right channel
     * @param random_dir Whether bouncing back should happen randomly
     * @param flat_gate Whether the gate is fully rectangular or the chambers are fully circular.
     */
    Simulation(int num_particles, double bridge_width, double circle_radius = 1.,
               double circle_distance = 0.5, int left_gate_capacity = 3,
               int right_gate_capacity = 3,
               bool random_dir = false, bool flat_gate = false);

    // Important parameters
    const int num_particles;
    int left_gate_capacity;
    int right_gate_capacity;
    /**
     * These values are used to compute the current.
     * The current is defined as the number of particles that move from left to right
     * (both through the center and the back channel) in a certain time interval, divided by that time.
     */
    std::vector<int> current_counters;
    const int FROM_LEFT_TO_RIGHT_INNER = 0;
    const int FROM_LEFT_TO_RIGHT_OUTER = 1;
    const int FROM_RIGHT_TO_LEFT_INNER = 2;
    const int FROM_RIGHT_TO_LEFT_OUTER = 3;
    unsigned long in_left;
    unsigned long num_collisions = 0;
    // Other parameters
    double circle_radius;
    double circle_distance;
    double bridge_width;
    double second_width;
    // Computed quantities
    double left_center_x;
    double right_center_x;
    double max_path;
    double bridge_length; // measured from the width
    double box_x_radius;
    double box_y_radius;
    // Choose second_length 0 to revert to old case
    double second_length; // also measured from the width;
    const unsigned long LEFT = 0;
    const unsigned long RIGHT = 1;
    bool explosion_direction_is_random;
    bool gate_is_flat;
    bool distance_as_channel_length = false;
    unsigned long expected_collisions = 0;

    // There is a (geometrical) difference between the distance between the urns and the length of the channel
    // if the gate is flat. While the former is nicer from a modelling point of view,
    // The latter provides an easier mathematical analysis.
    // To facilitate this, this boolean switch assumes circle distances as bridge lengths (and corrects for them)

    /**
     * Compute necessary parameters for the simulation, initialize data structures.
     * Run only once per simulation. Different runs require new setups and (therefore) new objects.
     */
    void setup();

    /**
     * Start the simulation. Initialize particles and times, and compute the next (first) impact.
     * @param left_ratio ratio of particles that should be initiated on the left side.
     */
    void start(double left_ratio);

    /**
     * Check if the point (x,y) is in the domain (collection for saying it is either in the left, the right or the gate.
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point
     * @return true if inside the domain, false otherwise.
     */
    bool is_in_domain(double x, double y) const;

    /**
     * Check if the point (x,y) is in a circle on side `side`.
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point
     * @param side `LEFT` or `RIGHT`
     * @return `true` if point in circle, `false` otherwise
     */
    bool is_in_circle(double x, double y, const unsigned long &side) const;

    /**
     * Check if point is in the bridge
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point
     * @return `true` if point in bridge, `false` otherwise
     */
    bool is_in_bridge(double x, double y) const;

    /**
     * Check if point is in the back channel
     * @param x x-coordinate of the point
     * @param y y-coordinate of the point
     * @return `true` if point in bridge, `false` otherwise
     */
    bool is_in_second_bridge(double x, double y) const;

    /**
     * Compute the intersections of a particle with an urn (on the x-axis)
     * @param particle Particle index
     * @param center_x x-coordinate of the urn
     * @param t1 First intersection
     * @param t2 Second intersection
     */
    void circle_intersections(const unsigned &particle, double center_x, double &t1, double &t2) const;

    /**
     * Compute the time it takes for a particle to reach the bridge.
     * If collision with bridge will not happen, return a time exceeding maximum collision time.
     * @param particle Particle index
     * @param normal_angle normal angle of the collision surface (output variable).
     * Only well-defined for not-exceeding-maximum-time collision time
     * @return time to next collision with bridge
     */
    double time_to_hit_bridge(const unsigned long &particle, double &normal_angle) const;

    /**
     * Compute the time it takes for a particle to reach the horizontal lines of the back channel.
     * If collision with bridge will not happen, return a time exceeding maximum collision time.
     * @param particle Particle index
     * @param normal_angle normal angle of the collision surface (output variable).
     * Only well-defined for not-exceeding-maximum-time collision time
     * @return time to next collision with bridge
     */
    double time_to_hit_second_bridge(const unsigned long &particle, double &normal_angle) const;

    /**
     * Computes the time it takes for a particle to reach the boundary of a reservoir
     * @param particle Particle index
     * @param center_x Center of the (circular) reservoir
     * @param normal_angle normal angle of the collision surface (output variable).
     * Only well-defined for not-exceeding-maximum-time collision time
     * @return time to next collision with circle
     */
    double time_to_hit_circle(const unsigned long &particle, double center_x, double &normal_angle) const;

    /**
     * Computes the time it takes for a particle to hit the gate.
     * No normal angle required because this collision does not reflect.
     * @param particle Particle index
     * @return time to next collision with gate
     */
    double time_to_hit_gate(const unsigned long &particle) const;

    /**
     * Time for a particle to hit the central vertical axis
     * @param particle Particle index
     * @return time to passing middle
     */
    double time_to_hit_middle(const unsigned long &particle) const;

    /**
     * Time for a particle to hit the outer vertical bounds
     * @param particle Particle index
     * @return time to passing middle
     */
    double time_to_hit_bounds(const unsigned long &particle) const;

    /**
     * Computes the next impact of a particle by finding the minimum impact time of all options.
     * @param particle Particle index
     */
    void compute_next_impact(const unsigned long &particle);

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
    std::vector<int> gate_capacities;

    std::vector<unsigned long> currently_in_left_gate;
    std::vector<unsigned long> currently_in_right_gate;
    std::vector<std::vector<unsigned long>> gate_contents;
    std::vector<std::vector<bool>> gate_arrays;

    /**
     * Compute the current position of a particle, based on the `time` variable
     * Interpolates on the assumption that the particle has no collision from its last update until this time
     * (which should be a valid assumption) and extrapolates linearly for an exact position.
     * @param particle Particle index
     * @param x output variable for x position
     * @param y output variable for y position
     */
    void get_current_position(const unsigned long &particle, double &x, double &y) const;


    /**
     * (Re)set the particle to some initial position. We also use this method if we lose a particle due to tricky
     * arithmetical errors (which happen once every 100000 instances)
     * @param particle Particle index
     * @param box_x_radius Horizontal radius of the box
     * @param box_y_radius Vertical radius of the box
     * @param direction LEFT or RIGHT
     */
    void reset_particle(const unsigned long &particle, const unsigned long &direction);

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

    /**
     * Check if a particle crossed the back channel and needs to have periodic boundary conditions applied.
     * @param particle particle index
     */
    void check_boundary_condition(const unsigned long &particle);

    /**
     * Check if a particle crossed the center vertical axis and needs to be counted as a switch.
     * @param particle Particle index
     */
    void count_first_gate_crossing(const unsigned long &particle);

    /**
     * Check if the particle is in the gate on side `direction`
     * @param x x-coordinate of particle position
     * @param y y-coordinate of particle position
     * @param direction Side of gate
     * @return True if particle in that side of the gate, false otherwise
     */
    bool is_in_gate(double x, double y, const unsigned long &direction) const;

    /**
     * Check if particle is moving towards the center or away from the center.
     * @param particle Particle index
     * @return True if particle moves towards the center, false otherwise
     */
    bool is_going_in(const unsigned long &particle) const;

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
    void update(double write_dt);

    /**
     * Print the current status of the simulation to stdout
     */
    void print_status() const;

    /**
     * Write the positions to file at a certain time, which *should* be between the current time and the next collision.
     * However, this is not enforced in this method.
     * @param time Time at which positions would be interpolated.
     */
    void write_positions_to_file(double time) const;

    /**
     * Write bounce_map
     */
    void write_bounce_map_to_file(const unsigned long &particle) const;

    /**
     * Compute the mass spread in the chamber.
     * The mass spread is defined as the number of particles in the right urn minus the number of particles in the left urn
     * divided by the total number of particles. This means that an qqual distribution of mass yields 0, while
     * the right urn full and the left empty yields a mass spread of one.
     *
     * @return mass spread as a double between -1 and 1.
     */
    double get_mass_spread() const;

    /**
     * Compute the reflection angle based on an ingoing angle and the normal angle of the surface.
     * @param angle_in incoming angle of the particle
     * @param normal_angle normal angle of the surface
     * @return outgoing angle of the particle
     */
    double get_reflection_angle(double angle_in, double normal_angle) const;

    /**
     * Compute the angle a particle should have after a gate explosion (reverse/random)
     * @param particle Particle index
     * @return angle of the particle after exploding the gate
     */
    double get_retraction_angle(const unsigned long &particle) const;

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

    /**
     * Sort the indices of the particles with respect to the next collision.
     * By using sorted indices, we speed op the simulation 4-5x.
     */
    void sort_indices();

    /**
     * Find the index of a particle in a sorted array.
     * This method uses binary search.
     * @param particle Particle index
     * @return Index of the particle where it fits.
     */
    unsigned long find_index(const unsigned long &particle) const;

    /**
     * Updates the position of a particle in the sorted index list, based on its new collision time.
     * @param particle Particle index
     * @param was_minimum If the particle was at the base of the list, no need to search its old position.
     */
    void reindex_particle(const unsigned long &particle, bool was_minimum);

    /**
     * Insert a particle in the list where it belongs, submethod.
     * @param particle
     */
    void insert_index(const unsigned long &particle);

    std::shared_ptr<std::random_device> rd;
    std::shared_ptr<std::mt19937> rng;
    std::shared_ptr<std::uniform_real_distribution<double>> unif_real;
    int reset_counter = 0;
    std::vector<unsigned long> sorted_indices;
};


#endif //TERRIER_SIMULATION_H
