#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>
#include "simulation.h"
#include <cmath>

BOOST_AUTO_TEST_SUITE(test_simulation)
    double eps = 1E-9;

    Simulation get_sim(int num_particles) {
        auto sim = Simulation(num_particles, 0.3);
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.gate_capacity = 1;
        return sim;
    }


    BOOST_AUTO_TEST_CASE(test_simulation_creation) {
        auto sim = get_sim(1000);
        BOOST_CHECK_EQUAL(sim.num_particles, 1000);
    }

    BOOST_AUTO_TEST_CASE(test_inside_methods) {
        auto sim = get_sim(100);
        sim.circle_radius = 1;
        sim.circle_distance = 0.5;
        sim.bridge_height = 0.1;
        sim.gate_capacity = 1;
        sim.setup();
        double x = 0;
        double y = 0;
        BOOST_CHECK(not sim.is_in_left_circle(x, y));
        BOOST_CHECK(not sim.is_in_right_circle(x, y));
        BOOST_CHECK(sim.is_in_bridge(x, y));
        BOOST_CHECK(sim.is_in_gate_radius(x, y));
        BOOST_CHECK(sim.is_in_domain(x, y));
        x = -2;
        y = 0.3;
        BOOST_CHECK(not sim.is_in_bridge(x, y));
        BOOST_CHECK(sim.is_in_left_circle(x, y));
        BOOST_CHECK(not sim.is_in_right_circle(x, y));
        BOOST_CHECK(sim.is_in_domain(x, y));
        BOOST_CHECK(not sim.is_in_gate_radius(x, y));
        x = 2.25;
        y = 0;
        BOOST_CHECK(not sim.is_in_left_circle(x, y));
        y = 0.06;
        BOOST_CHECK(not(sim.is_in_domain(x, y) and sim.is_in_gate_radius(x, y)));
        x = 0;
        BOOST_CHECK(not sim.is_in_domain(x, y));
        BOOST_CHECK(sim.is_in_gate_radius(x, y));
    }

    BOOST_AUTO_TEST_CASE(test_particle_init) {
        auto sim = get_sim(1000);
        sim.setup();
        sim.start();
        BOOST_CHECK_EQUAL(sim.total_left.at(0), 1000);
        BOOST_CHECK_EQUAL(sim.total_right.at(0), 0);
        bool correct = true;
        for (int i = 0; i < 1000; i++) {
            correct &= sim.is_in_left_circle(sim.positions(i, 0), sim.positions(i, 1));
        }
        BOOST_CHECK(correct);
    }


    BOOST_AUTO_TEST_CASE(test_bridge_collision) {
        auto sim = get_sim(1);
        double sqrt_2 = std::sqrt(2);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();

        //When lines intersect, we want a hit
        sim.positions(0, 0) = 0;
        sim.positions(0, 1) = 0;
        sim.directions(0) = pi / 2; // straight up
        double angle = 2;
        double time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_height / 2, eps);
        BOOST_CHECK_CLOSE(angle, -pi / 2, eps);
        sim.directions(0) = pi / 4; // right up
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_height * sqrt_2 / 2, eps);
        sim.positions(0, 1) = -.1;
        sim.directions(0) = pi * 3 / 4; // left up
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_height / 2 * sqrt_2, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, eps);

        // When no intersection, we want no hit
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = -pi / 2; // straight up
        angle = 2;
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        sim.positions(0, 0) = sim.right_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = -pi * 3 / 4; // straight up
        angle = 2;
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        // Touching is also intersecting
        sim.positions(0, 0) = -sim.bridge_size / 2 - 0.1;
        sim.positions(0, 1) = sim.bridge_height / 2 + 0.1;
        sim.directions(0) = -pi * 1 / 4; // from top left to touch top left bridge point
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, 0.1 * sqrt_2, eps);
        BOOST_CHECK_CLOSE(angle, -pi / 2, eps);
        sim.positions(0, 0) = -sim.bridge_size / 2 - 0.1;
        sim.positions(0, 1) = sim.bridge_height / 2 + 0.1;
        sim.directions(0) = 2 * pi; // Straight right
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        //Parallel is never intersecting
        sim.positions(0, 1) = 0;
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        // Not even if the line segments overlap
        sim.positions(0, 0) = -sim.bridge_size;
        sim.positions(0, 1) = sim.bridge_height / 2;
        sim.directions(0) = 2 * pi; // Straight right
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        // What if we leave from the bridge? Then we don't want a hit
        sim.positions(0, 0) = sim.bridge_size / 3;
        sim.positions(0, 1) = sim.bridge_height / 2;
        sim.directions(0) = -pi / 2; // Straight down
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_height, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, eps);
    }

    BOOST_AUTO_TEST_CASE(test_circle_collision) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        double time, angle;
        // standard: from center of circle to boundary
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = -pi / 2;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, eps);
        sim.directions(0) = pi / 4;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle, -pi * 3 / 4, eps);
        sim.directions(0) = -pi;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle + 2 * pi, 0 + 2 * pi, eps); // weird test...
        sim.positions(0, 0) = sim.right_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = -pi / 2;
        time = sim.time_to_hit_circle(0, sim.right_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, 0.01);
        //Inscribe the circle with a equilateral triangle, top up.
        // Now go from the left down point to the other two.
        double x = 1.5 / sqrt(3) * sim.circle_radius;
        sim.positions(0, 0) = sim.left_center_x - x;
        sim.positions(0, 1) = -sim.circle_radius / 2;
        sim.directions(0) = pi / 3;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, 2 * x, eps);
        BOOST_CHECK_CLOSE(angle, -pi / 2, eps);
        sim.directions(0) = 2 * pi;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, 2 * x, eps);
        BOOST_CHECK_CLOSE(angle, pi * 5 / 6, eps);
    }

    BOOST_AUTO_TEST_CASE(test_circle_consistency) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        sim.update();
        BOOST_CHECK(sim.is_in_domain(sim.positions(0, 0), sim.positions(0, 1)));
        sim.update();
        BOOST_CHECK(sim.is_in_domain(sim.positions(0, 0), sim.positions(0, 1)));
        sim.update();
        BOOST_CHECK(sim.is_in_domain(sim.positions(0, 0), sim.positions(0, 1)));
    }

    BOOST_AUTO_TEST_CASE(test_bridge_coupling) {
        auto sim = get_sim(1);
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        BOOST_CHECK(sim.bridge_size > sim.circle_distance);
        BOOST_CHECK(sim.is_in_domain(sim.bridge_size / 2 - 0.001, sim.bridge_height / 2 - 0.001));
        BOOST_CHECK(not sim.is_in_right_circle(sim.bridge_size / 2 - 0.001, sim.bridge_height / 2 - 0.001));
    }

    BOOST_AUTO_TEST_CASE(test_circle_bridge_connection) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        double time, angle;
        // standard: from center of circle to boundary
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = 2 * pi;
        time = sim.time_to_hit_circle(0, sim.right_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius * 3 + 0.5, eps);
        BOOST_CHECK_CLOSE(angle, pi, eps);
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);
    }

    BOOST_AUTO_TEST_CASE(test_hit_gate) {
        auto sim = get_sim(1);
        sim.gate_radius = 0.6;
        double pi = 3.141592653589793;
        sim.setup();
        sim.start();
        double time, angle;
        // standard: from center of circle to boundary
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = 2 * pi;
        time = sim.time_to_hit_gate(0);
        BOOST_CHECK_CLOSE(time, -sim.left_center_x - sim.gate_radius, eps);
        // From below
        sim.positions(0, 0) = 0;
        sim.positions(0, 1) = -4;
        sim.directions(0) = pi / 2;
        time = sim.time_to_hit_gate(0);
        BOOST_CHECK_CLOSE(time, 4 - sim.gate_radius, eps);
        // From inside
        sim.positions(0, 0) = -sim.gate_radius;
        sim.positions(0, 1) = 0;
        sim.directions(0) = 2 * pi;
        time = sim.time_to_hit_gate(0);
        BOOST_CHECK_CLOSE(time, 2 * sim.gate_radius, eps);
    }

    BOOST_AUTO_TEST_CASE(test_updates) {
        auto sim = get_sim(10000);
        sim.setup();
        sim.start();
        sim.measure();
        BOOST_CHECK(sim.in_left.sum() + sim.in_right.sum() <= sim.num_particles);
    }

    BOOST_AUTO_TEST_CASE(test_reflection_angle) {
        auto sim = get_sim(10000);
        double pi = 3.141592653589793;
        // pretty self-explanatory. Do the math
        double angle_in = pi / 6;
        double normal_angle = pi;
        BOOST_CHECK_CLOSE(sim.get_reflection_angle(angle_in, normal_angle), pi * 5 / 6, eps);
        angle_in = pi / 2;
        normal_angle = -pi / 2;
        BOOST_CHECK_CLOSE(sim.get_reflection_angle(angle_in, normal_angle), -pi / 2, eps);
        angle_in = 0;
        normal_angle = pi * 5 / 4;
        BOOST_CHECK_CLOSE(sim.get_reflection_angle(angle_in, normal_angle), 3 * pi / 2, eps);
        angle_in = 0;
        normal_angle = pi / 2;
        BOOST_CHECK_CLOSE(sim.get_reflection_angle(angle_in, normal_angle), 0, eps);
    }

    BOOST_AUTO_TEST_CASE(test_collision) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        // standard: from center of circle to bottom boundary
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = -pi / 2;
        sim.compute_next_impact(0);
        BOOST_CHECK_CLOSE(sim.next_impact_times(0), sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(sim.next_directions(0), pi / 2, eps);
        // from center of circle to
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.directions(0) = 0;
        sim.gate_radius = 0; // don't do in real life
        sim.compute_next_impact(0);
        BOOST_CHECK_CLOSE(sim.next_impact_times(0), sim.circle_radius * 3 + 0.5, eps);
        BOOST_CHECK_CLOSE(sim.next_directions(0), pi, eps);
    }

    BOOST_AUTO_TEST_CASE(test_update) {
        // One of the most important tests, I reckon
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        // standard: from center of circle to bottom boundary
        sim.positions(0, 0) = sim.left_center_x;
        sim.positions(0, 1) = 0;
        sim.next_directions(0) = -pi / 2;
        sim.next_impact_times(0) = 0;
        sim.next_positions(0, 0) = sim.left_center_x;
        sim.next_positions(0, 1) = 0;

        sim.update();
        BOOST_CHECK_CLOSE(sim.next_impact_times(0), sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(sim.next_directions(0), pi / 2, eps);
        sim.update();
        BOOST_CHECK_CLOSE(sim.positions(0, 1), -sim.circle_radius, eps);
    }

    BOOST_AUTO_TEST_CASE(test_collective_update) {
        // One of the most important tests, I reckon
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        // Make an inscribed square, in the direction of the clock
        double side = sim.circle_radius / sqrt(2);
        sim.positions(0, 0) = -sim.left_center_x - side;
        sim.positions(0, 1) = 0;
        sim.directions(0) = pi / 2;
        sim.compute_next_impact(0);
        sim.update();
        BOOST_CHECK_CLOSE(sim.positions(0, 0), -sim.left_center_x - side, eps);
        BOOST_CHECK_CLOSE(sim.positions(0, 1), side, eps);
        BOOST_CHECK_CLOSE(sim.directions(0), 0, eps);
        sim.update();
        BOOST_CHECK_CLOSE(sim.positions(0, 0), -sim.left_center_x + side, eps);
        BOOST_CHECK_CLOSE(sim.positions(0, 1), side, eps);
        BOOST_CHECK_CLOSE(sim.directions(0), -pi / 2, eps);
        sim.update();
        BOOST_CHECK_CLOSE(sim.positions(0, 0), -sim.left_center_x + side, eps);
        BOOST_CHECK_CLOSE(sim.positions(0, 1), -side, eps);
        BOOST_CHECK_CLOSE(sim.directions(0), pi, eps);
        sim.update();
        BOOST_CHECK_CLOSE(sim.positions(0, 0), -sim.left_center_x - side, eps);
        BOOST_CHECK_CLOSE(sim.positions(0, 1), -side, eps);
        BOOST_CHECK_CLOSE(sim.directions(0), pi / 2, eps);
    }


    BOOST_AUTO_TEST_CASE(test_long_term_consistency) {
        auto sim = get_sim(1000);
        double pi = 3.141592653589793;
        sim.bridge_height = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start();
        bool all_correct = true;
        while (sim.time < 100) {
            sim.update();
            for (int i = 0; i < sim.num_particles; i++) {
                all_correct = all_correct and sim.is_in_domain(sim.positions(i, 0), sim.positions(i, 1));
            }
            BOOST_CHECK(all_correct);
        }
    }

    BOOST_AUTO_TEST_CASE(test_retraction_angle) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.circle_distance = 0.5;
        sim.gate_radius = 0.3;
        sim.setup();
        sim.start();
        // in left
        sim.positions(0, 0) = -0.27;
        sim.positions(0, 1) = 0.02;
        sim.directions(0) = -0.1;
        double angle = sim.get_retraction_angle(0);
        BOOST_CHECK(angle > pi / 2 and angle < pi * 3 / 2);
        // in right
        sim.positions(0, 0) = 0.26;
        sim.directions(0) = -pi;
        angle = sim.get_retraction_angle(0);
        BOOST_CHECK(abs(angle) < pi);
        // in bridge
        sim.positions(0, 0) = -0.22;
        sim.directions(0) = -pi;
        angle = sim.get_retraction_angle(0);
        BOOST_CHECK(abs(angle) < pi);
        // past bridge
        sim.positions(0, 0) = 0.22;
        sim.positions(0, 1) = 0.22;
        sim.directions(0) = pi / 6;
        angle = sim.get_retraction_angle(0);
        BOOST_CHECK(cos(angle) < 0);
        // or, more precisely
        BOOST_CHECK_CLOSE(angle, -pi * 3 / 4, eps);

    }




BOOST_AUTO_TEST_SUITE_END();