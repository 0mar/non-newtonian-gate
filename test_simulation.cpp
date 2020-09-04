#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>
#include "simulation.h"
#include <cmath>

BOOST_AUTO_TEST_SUITE(test_simulation)
    double eps = 1E-9;

    Simulation get_sim(int num_particles) {
        auto sim = Simulation(num_particles, 0.1);
        sim.circle_distance = 0.5;
        sim.left_gate_capacity = 1;
        sim.right_gate_capacity = 1;
        sim.explosion_direction_is_random = false;
        sim.gate_is_flat = true;
        return sim;
    }


    BOOST_AUTO_TEST_CASE(test_simulation_creation) {
        auto sim = get_sim(1000);
        BOOST_CHECK_EQUAL(sim.num_particles, 1000);
    }

    BOOST_AUTO_TEST_CASE(test_mass_spread) {
        auto sim = get_sim(1000);
        sim.setup();
        sim.start(0.5);
        BOOST_CHECK_CLOSE(sim.get_mass_spread(), 0, eps);
        auto sim2 = get_sim(100);
        sim2.setup();
        sim2.start(.25);
        BOOST_CHECK_CLOSE(std::fabs(sim2.get_mass_spread()), .5, eps);
        auto sim3 = get_sim(50);
        sim3.setup();
        sim3.start(1);
        BOOST_CHECK_CLOSE(std::fabs(sim3.get_mass_spread()), 1, eps);
    }

    BOOST_AUTO_TEST_CASE(test_bridge_length) {
        auto sim = get_sim(1);
        sim.circle_radius = 1;
        sim.circle_distance = 0.5;
        sim.bridge_width = 0.1;
        sim.setup();
        const double cutoff_radius = std::sqrt(std::pow(sim.circle_radius, 2) - std::pow(sim.bridge_width, 2) / 4);
        BOOST_CHECK_CLOSE(sim.circle_distance, sim.bridge_length + 2 * cutoff_radius - 2 * sim.circle_radius, eps);
    }

    BOOST_AUTO_TEST_CASE(test_new_bridge_coupling) {
        auto sim = get_sim(1);
        sim.distance_as_channel_length = false;
        sim.setup();
        double lx = 0;
        double ly = sim.bridge_width / 2 + eps * 0.01;
        double angle;
        sim.x_pos.at(0) = lx;
        sim.y_pos.at(0) = ly;
        sim.directions.at(0) = 0;
        BOOST_CHECK_CLOSE(sim.bridge_length, 2 * sim.time_to_hit_circle(0, sim.right_center_x, angle), eps);
    }

    BOOST_AUTO_TEST_CASE(test_bridge_compensation) {
        auto sim = get_sim(1);
        const double old_circle_distance = sim.circle_distance;
        sim.distance_as_channel_length = true;
        sim.setup();
        BOOST_CHECK_CLOSE(old_circle_distance, sim.bridge_length, eps);
        BOOST_CHECK(sim.circle_distance < sim.bridge_length);
    }

    BOOST_AUTO_TEST_CASE(test_bridge_coupling) {
        auto sim = get_sim(1);
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        BOOST_CHECK(sim.bridge_length > sim.circle_distance);
        BOOST_CHECK(sim.is_in_domain(sim.bridge_length / 2 - 0.001, sim.bridge_width / 2 - 0.001));
        BOOST_CHECK(not sim.is_in_circle(sim.bridge_length / 2 - 0.001, sim.bridge_width / 2 - 0.001, sim.RIGHT));
    }

    BOOST_AUTO_TEST_CASE(test_inside_methods) {
        auto sim = get_sim(100);
        sim.circle_radius = 1;
        sim.circle_distance = 0.5;
        sim.bridge_width = 0.1;
        sim.left_gate_capacity = 1;
        sim.setup();
        double x = 0;
        double y = 0;
        BOOST_CHECK(not sim.is_in_circle(x, y, sim.LEFT));
        BOOST_CHECK(not sim.is_in_circle(x, y, sim.RIGHT));
        BOOST_CHECK(sim.is_in_bridge(x, y));
        BOOST_CHECK(sim.is_in_gate(x, y, sim.LEFT));
        BOOST_CHECK(sim.is_in_gate(x, y, sim.RIGHT));
        BOOST_CHECK(sim.is_in_domain(x, y));
        x = -2;
        y = 0.3;
        BOOST_CHECK(not sim.is_in_bridge(x, y));
        BOOST_CHECK(sim.is_in_circle(x, y, sim.LEFT));
        BOOST_CHECK(not sim.is_in_circle(x, y, sim.RIGHT));
        BOOST_CHECK(sim.is_in_domain(x, y));
        BOOST_CHECK(not sim.is_in_gate(x, y, sim.LEFT));
        BOOST_CHECK(not sim.is_in_gate(x, y, sim.RIGHT));
        x = 2.25;
        y = 0;
        BOOST_CHECK(not sim.is_in_circle(x, y, sim.LEFT));
        y = 0.06;
        BOOST_CHECK(not(sim.is_in_domain(x, y) and sim.is_in_gate(x, y, sim.RIGHT)));
        x = 0;
        BOOST_CHECK(not sim.is_in_domain(x, y));
        BOOST_CHECK(sim.is_in_gate(x, y, sim.LEFT));
        BOOST_CHECK(sim.is_in_gate(x, y, sim.RIGHT));
    }

    BOOST_AUTO_TEST_CASE(test_particle_init) {
        auto sim = get_sim(1000);
        sim.setup();
        sim.start(1);
        BOOST_CHECK_EQUAL(sim.num_particles, 1000);
        bool correct = true;
        for (unsigned long i = 0; i < 1000; i++) {
            correct &= sim.is_in_circle(sim.x_pos.at(i), sim.y_pos.at(i), sim.LEFT);
        }
        BOOST_CHECK(correct);
    }

    BOOST_AUTO_TEST_CASE(test_in_circle) {
        auto sim = get_sim(1);
        double radius = 0.5;
        sim.circle_distance = radius * 2;
        sim.setup();
        sim.start(1);
        sim.gate_is_flat = false;
        // test the flat and nonflat cases
        for (unsigned int i = 0; i < 2; i++) {
            if (i == 1) {
                radius = sim.bridge_length / 2;
                sim.gate_is_flat = true;
            }
            BOOST_CHECK(sim.is_in_gate(-radius + eps, 0, sim.LEFT));
            BOOST_CHECK(not sim.is_in_gate(-radius - eps, 0, sim.RIGHT));
            BOOST_CHECK(sim.is_in_gate(+radius - eps, 0, sim.RIGHT));
            BOOST_CHECK(not sim.is_in_gate(+radius + eps, 0, sim.LEFT));
            BOOST_CHECK(sim.is_in_gate(-radius / 2, 0.05, sim.LEFT));
        }
    }

    BOOST_AUTO_TEST_CASE(test_to_middle) {
        auto sim = get_sim(1);
        sim.setup();
        sim.start(1);
        sim.x_pos.at(0) = -.5;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = 0;
        double distance = sim.time_to_hit_middle(0);
        BOOST_CHECK_CLOSE(distance, 0.5, eps);
    }

    BOOST_AUTO_TEST_CASE(test_bridge_collision) {
        auto sim = get_sim(1);
        double sqrt_2 = std::sqrt(2);
        double pi = 3.141592653589793;
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);

        //When lines intersect, we want a hit
        sim.x_pos.at(0) = 0;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = pi / 2; // straight up
        double angle = 2;
        double time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_width / 2, eps);
        BOOST_CHECK_CLOSE(angle, -pi / 2, eps);
        sim.directions.at(0) = pi / 4; // right up
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_width * sqrt_2 / 2, eps);
        sim.y_pos.at(0) = -.1;
        sim.directions.at(0) = pi * 3 / 4; // left up
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_width / 2 * sqrt_2, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, eps);

        // When no intersection, we want no hit
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = -pi / 2; // straight up
        angle = 2;
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        sim.x_pos.at(0) = sim.right_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = -pi * 3 / 4; // straight up
        angle = 2;
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        // Touching is also intersecting
        sim.x_pos.at(0) = -sim.bridge_length / 2 - 0.1;
        sim.y_pos.at(0) = sim.bridge_width / 2 + 0.1;
        sim.directions.at(0) = -pi * 1 / 4; // from top left to touch top left bridge point
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, 0.1 * sqrt_2, eps);
        BOOST_CHECK_CLOSE(angle, -pi / 2, eps);
        sim.x_pos.at(0) = -sim.bridge_length / 2 - 0.1;
        sim.y_pos.at(0) = sim.bridge_width / 2 + 0.1;
        sim.directions.at(0) = 2 * pi; // Straight right
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        //Parallel is never intersecting
        sim.y_pos.at(0) = 0;
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        // Not even if the line segments overlap
        sim.x_pos.at(0) = -sim.bridge_length;
        sim.y_pos.at(0) = sim.bridge_width / 2;
        sim.directions.at(0) = 2 * pi; // Straight right
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);

        // What if we leave from the bridge? Then we don't want a hit
        sim.x_pos.at(0) = sim.bridge_length / 3;
        sim.y_pos.at(0) = sim.bridge_width / 2;
        sim.directions.at(0) = -pi / 2; // Straight down
        time = sim.time_to_hit_bridge(0, angle);
        BOOST_CHECK_CLOSE(time, sim.bridge_width, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, eps);
    }

    BOOST_AUTO_TEST_CASE(test_circle_collision) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        double time, angle;
        // standard: from center of circle to boundary
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = -pi / 2;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, eps);
        sim.directions.at(0) = pi / 4;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle, -pi * 3 / 4, eps);
        sim.directions.at(0) = -pi;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle + 2 * pi, 0 + 2 * pi, eps); // weird test...
        sim.x_pos.at(0) = sim.right_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = -pi / 2;
        time = sim.time_to_hit_circle(0, sim.right_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(angle, pi / 2, 0.01);
        //Inscribe the circle with a equilateral triangle, top up.
        // Now go from the left down point to the other two.
        double x = 1.5 / sqrt(3) * sim.circle_radius;
        sim.x_pos.at(0) = sim.left_center_x - x;
        sim.y_pos.at(0) = -sim.circle_radius / 2;
        sim.directions.at(0) = pi / 3;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, 2 * x, eps);
        BOOST_CHECK_CLOSE(angle, -pi / 2, eps);
        sim.directions.at(0) = 2 * pi;
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, 2 * x, eps);
        BOOST_CHECK_CLOSE(angle, pi * 5 / 6, eps);
    }

    BOOST_AUTO_TEST_CASE(test_circle_consistency) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        sim.update(0);
        BOOST_CHECK(sim.is_in_domain(sim.x_pos.at(0), sim.y_pos.at(0)));
        sim.update(0);
        BOOST_CHECK(sim.is_in_domain(sim.x_pos.at(0), sim.y_pos.at(0)));
        sim.update(0);
        BOOST_CHECK(sim.is_in_domain(sim.x_pos.at(0), sim.y_pos.at(0)));
    }

    BOOST_AUTO_TEST_CASE(test_circle_bridge_connection) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        double time, angle;
        // standard: from center of circle to boundary
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = 2 * pi;
        time = sim.time_to_hit_circle(0, sim.right_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.circle_radius * 3 + 0.5, eps);
        BOOST_CHECK_CLOSE(angle, pi, eps);
        time = sim.time_to_hit_circle(0, sim.left_center_x, angle);
        BOOST_CHECK_CLOSE(time, sim.max_path, eps);
    }

    BOOST_AUTO_TEST_CASE(test_hit_gate) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.setup();
        sim.start(1);
        double time, angle;
        // standard: from center of circle to boundary
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = 2 * pi;
        sim.gate_is_flat = true;
        time = sim.time_to_hit_gate(0);

        BOOST_CHECK_CLOSE(time, -sim.left_center_x - sim.bridge_length / 2, eps);
        sim.gate_is_flat = false;
        time = sim.time_to_hit_gate(0);
        BOOST_CHECK_CLOSE(time, -sim.left_center_x - sim.circle_distance / 2, eps);

        // From inside
        sim.gate_is_flat = true;
        sim.x_pos.at(0) = -sim.bridge_length / 2;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = 2 * pi;
        time = sim.time_to_hit_gate(0);
        BOOST_CHECK_CLOSE(time, sim.bridge_length, eps);
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
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        // standard: from center of circle to bottom boundary
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = -pi / 2;
        sim.compute_next_impact(0);
        sim.compute_next_impact(0);
        BOOST_CHECK_CLOSE(sim.next_impact_times.at(0), sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(sim.next_directions.at(0), pi / 2, eps);
        // from center of circle to
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = 0;
//        sim.gate_radius = 0; // don't do in real life
//        sim.compute_next_impact(0);
//        sim.compute_next_impact(0);
//        BOOST_CHECK_CLOSE(sim.next_impact_times.at(0), sim.circle_radius * 3 + 0.5, eps);
//        BOOST_CHECK_CLOSE(sim.next_directions.at(0), pi, eps);
    }

    BOOST_AUTO_TEST_CASE(test_update) {
        // One of the most important tests, I reckon
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        // standard: from center of circle to bottom boundary
        sim.x_pos.at(0) = sim.left_center_x;
        sim.y_pos.at(0) = 0;
        sim.next_directions.at(0) = -pi / 2;
        sim.next_impact_times.at(0) = 0;
        sim.next_x_pos.at(0) = sim.left_center_x;
        sim.next_y_pos.at(0) = 0;

        sim.update(0);
        BOOST_CHECK_CLOSE(sim.next_impact_times.at(0), sim.circle_radius, eps);
        BOOST_CHECK_CLOSE(sim.next_directions.at(0), pi / 2, eps);
        sim.update(0);
        BOOST_CHECK_CLOSE(sim.y_pos.at(0), -sim.circle_radius, eps);
    }

    BOOST_AUTO_TEST_CASE(test_collective_update) {
        // One of the most important tests, I reckon
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        // Make an inscribed square, in the direction of the clock
        double side = sim.circle_radius / sqrt(2);
        sim.x_pos.at(0) = -sim.left_center_x - side;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = pi / 2;
        sim.compute_next_impact(0);
        sim.update(0);
        BOOST_CHECK_CLOSE(sim.x_pos.at(0), -sim.left_center_x - side, eps);
        BOOST_CHECK_CLOSE(sim.y_pos.at(0), side, eps);
        BOOST_CHECK_CLOSE(sim.directions.at(0), 0, eps);
        sim.update(0);
        BOOST_CHECK_CLOSE(sim.x_pos.at(0), -sim.left_center_x + side, eps);
        BOOST_CHECK_CLOSE(sim.y_pos.at(0), side, eps);
        BOOST_CHECK_CLOSE(sim.directions.at(0), -pi / 2, eps);
        sim.update(0);
        BOOST_CHECK_CLOSE(sim.x_pos.at(0), -sim.left_center_x + side, eps);
        BOOST_CHECK_CLOSE(sim.y_pos.at(0), -side, eps);
        BOOST_CHECK_CLOSE(sim.directions.at(0), pi, eps);
        sim.update(0);
        BOOST_CHECK_CLOSE(sim.x_pos.at(0), -sim.left_center_x - side, eps);
        BOOST_CHECK_CLOSE(sim.y_pos.at(0), -side, eps);
        BOOST_CHECK_CLOSE(sim.directions.at(0), pi / 2, eps);
    }


    BOOST_AUTO_TEST_CASE(test_long_term_consistency) {
        // Takes a long time
        auto sim = get_sim(1000);
        double pi = 3.141592653589793;
        sim.bridge_width = 0.1;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        bool all_correct = true;
        bool left_or_right = true;
        while (sim.time < 40) {
            sim.update(0);
            for (unsigned long i = 0; i < sim.num_particles; i++) {
                all_correct = all_correct and sim.is_in_domain(sim.x_pos.at(i), sim.y_pos.at(i));
            }
            BOOST_CHECK(all_correct);
            BOOST_CHECK(left_or_right);
        }
    }

    BOOST_AUTO_TEST_CASE(test_retraction_angle) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.circle_distance = 0.5;
        sim.setup();
        sim.start(1);
        // in left
        sim.x_pos.at(0) = -0.27;
        sim.y_pos.at(0) = 0.02;
        sim.directions.at(0) = -0.1;
        double angle = sim.get_retraction_angle(0);
        if (sim.explosion_direction_is_random) {
            BOOST_CHECK(angle > pi / 2 and angle < pi * 3 / 2);
        } else {
            BOOST_CHECK_CLOSE(angle, 0.1 + pi, eps);
        }

        // in right
        sim.x_pos.at(0) = 0.26;
        sim.directions.at(0) = -pi;
        angle = sim.get_retraction_angle(0);
        if (sim.explosion_direction_is_random) {
            BOOST_CHECK(std::abs(angle) < pi);
        } else {
            BOOST_CHECK_CLOSE(angle - pi, sim.directions.at(0) + 2 * pi, eps);
        }
        // in bridge
        sim.x_pos.at(0) = -0.22;
        sim.directions.at(0) = -pi;
        angle = sim.get_retraction_angle(0);
        if (sim.explosion_direction_is_random) {
            BOOST_CHECK(angle > pi / 2 and angle < pi * 3 / 2);
        } else {
            BOOST_CHECK_CLOSE(angle, sim.directions.at(0), eps);
        }
        // past bridge
        sim.x_pos.at(0) = 0.22;
        sim.y_pos.at(0) = 0.22;
        sim.directions.at(0) = pi / 6;
        angle = sim.get_retraction_angle(0);
        if (sim.explosion_direction_is_random) {
            BOOST_CHECK(cos(angle) > 0);
        } else {
            BOOST_CHECK_CLOSE(angle, sim.directions.at(0), eps);
        }
    }

    BOOST_AUTO_TEST_CASE(test_get_current_position) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.setup();
        sim.start(1);
        sim.x_pos.at(0) = 0;
        sim.y_pos.at(0) = 0;
        sim.time = sim.next_impact_times.at(0) / 2;
        double x, y;
        sim.get_current_position(0, x, y);
        BOOST_CHECK_CLOSE(x, sim.next_x_pos.at(0) / 2, eps);
        BOOST_CHECK_CLOSE(y, sim.next_y_pos.at(0) / 2, eps);

    }

    BOOST_AUTO_TEST_CASE(test_going_in_out_flat_gate) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.circle_radius = 1;
        sim.circle_distance = 0.5;
        sim.bridge_width = 0.3;
        sim.gate_is_flat = true;
        // Hitting from a position radially to the gate.
        // Pick some angle small enough
        double angle = 1. / 12 * pi;
        sim.setup();
        sim.start(1);
        sim.x_pos.at(0) = sim.circle_radius * cos(pi - angle);
        sim.y_pos.at(0) = sim.circle_radius * sin(pi - angle);
        BOOST_CHECK(sim.is_in_circle(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(sim.gate_contents.at(sim.LEFT).empty());
        sim.directions.at(0) = -angle;
        sim.compute_next_impact(0);
        double next_x = -sim.bridge_length / 2;
        double next_y = sim.bridge_length / 2 * tan(angle);
        sim.update(0);
        BOOST_CHECK_CLOSE(next_x, sim.x_pos.at(0), eps);
        BOOST_CHECK_CLOSE(next_y, sim.y_pos.at(0), eps);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT) and
                    not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.RIGHT));
        BOOST_CHECK_CLOSE(sim.directions.at(0), -angle, eps);
        BOOST_CHECK(sim.gate_contents.at(sim.LEFT).size() == 1);
        BOOST_CHECK(sim.gate_contents.at(sim.RIGHT).empty());
        sim.update(0);
        BOOST_CHECK_CLOSE(eps, pow(sim.x_pos.at(0), 2) + eps, eps); // more weird unclose
        BOOST_CHECK_CLOSE(eps, pow(sim.y_pos.at(0), 2) + eps, eps);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.RIGHT) and
                    not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK_CLOSE(sim.directions.at(0), -angle, eps);
        BOOST_CHECK(sim.gate_contents.at(sim.RIGHT).empty()); // Discount leaving particles
        BOOST_CHECK(sim.gate_contents.at(sim.LEFT).empty());
        sim.update(0);
        BOOST_CHECK_CLOSE(-next_x, sim.x_pos.at(0), eps);
        BOOST_CHECK_CLOSE(-next_y, sim.y_pos.at(0), eps);
        BOOST_CHECK(not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.RIGHT));
        BOOST_CHECK_CLOSE(sim.directions.at(0), -angle, eps);
        BOOST_CHECK_EQUAL(sim.gate_contents.at(sim.LEFT).size() + sim.gate_contents.at(sim.RIGHT).size(), 0);

    }

    BOOST_AUTO_TEST_CASE(test_going_in_out_hollow_gate) {
        auto sim = get_sim(1);
        double pi = 3.141592653589793;
        sim.circle_radius = 1;
        sim.circle_distance = 0.5;
        sim.bridge_width = 0.3;
        sim.gate_is_flat = false;
        // Hitting from a position radially to the gate.
        // Pick some angle small enough
        double angle = 1. / 12 * pi;
        sim.setup();
        sim.start(1);
        sim.x_pos.at(0) = sim.circle_radius * cos(pi - angle);
        sim.y_pos.at(0) = sim.circle_radius * sin(pi - angle);
        BOOST_CHECK(sim.is_in_circle(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(sim.gate_contents.at(sim.LEFT).empty());
        sim.directions.at(0) = -angle;
        sim.compute_next_impact(0);
        sim.update(0);
        double next_x = sim.x_pos.at(0);
        double next_y = sim.y_pos.at(0);
        double radius_sqr = (next_x - sim.left_center_x) * (next_x - sim.left_center_x) + next_y * next_y;
        BOOST_CHECK_CLOSE(radius_sqr, sim.circle_radius * sim.circle_radius, eps);
        BOOST_CHECK(next_y > 0);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT) and
                    not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.RIGHT));
        BOOST_CHECK_CLOSE(sim.directions.at(0), -angle, eps);
        BOOST_CHECK(sim.gate_contents.at(sim.LEFT).size() == 1);
        BOOST_CHECK(sim.gate_contents.at(sim.RIGHT).empty());
        sim.update(0);
        BOOST_CHECK_CLOSE(eps, pow(sim.x_pos.at(0), 2) + eps, eps); // more weird unclose
        BOOST_CHECK_CLOSE(eps, pow(sim.y_pos.at(0), 2) + eps, eps);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.RIGHT) and
                    not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK_CLOSE(sim.directions.at(0), -angle, eps);
        BOOST_CHECK(sim.gate_contents.at(sim.RIGHT).empty()); // Discount leaving particles
        BOOST_CHECK(sim.gate_contents.at(sim.LEFT).empty());
        sim.update(0);
        next_x = sim.x_pos.at(0);
        next_y = sim.y_pos.at(0);
        radius_sqr = (next_x - sim.right_center_x) * (next_x - sim.right_center_x) + next_y * next_y;
        BOOST_CHECK_CLOSE(radius_sqr, sim.circle_radius * sim.circle_radius, eps);
        BOOST_CHECK(next_y < 0);
        BOOST_CHECK(not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.RIGHT));
        BOOST_CHECK_CLOSE(sim.directions.at(0), -angle, eps);
        BOOST_CHECK_EQUAL(sim.gate_contents.at(sim.LEFT).size() + sim.gate_contents.at(sim.RIGHT).size(), 0);

    }


    BOOST_AUTO_TEST_CASE(test_is_going_in) {
        auto sim = get_sim(1);
        sim.setup();
        sim.start(1);
        sim.x_pos.at(0) = -sim.circle_radius;
        sim.y_pos.at(0) = 0;
        double angle = 0.1;
        const double pi = 3.141592653589793;
        sim.directions.at(0) = angle;
        BOOST_CHECK(sim.is_going_in(0));
        sim.directions.at(0) = -angle;
        BOOST_CHECK(sim.is_going_in(0));
        sim.directions.at(0) = pi - angle;
        BOOST_CHECK(not sim.is_going_in(0));
        sim.x_pos.at(0) = sim.circle_radius;
        BOOST_CHECK(sim.is_going_in(0));
    }

    BOOST_AUTO_TEST_CASE(test_explosion_bug) {
        auto sim = get_sim(1);
        sim.gate_is_flat = true;
        sim.setup();
        sim.start(1);
        double delta = 0.001;
        const double pi = 3.141592653589793;
        sim.x_pos.at(0) = -sim.bridge_length / 2 - delta * delta;
        sim.y_pos.at(0) = 0;
        sim.directions.at(0) = pi / 2 - delta;
        BOOST_CHECK(sim.is_going_in(0));
        BOOST_CHECK(not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(not sim.is_in_bridge(sim.x_pos.at(0), sim.y_pos.at(0)));
        sim.compute_next_impact(0);
        sim.update(0);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(sim.is_in_bridge(sim.x_pos.at(0), sim.y_pos.at(0)));
        sim.update(0);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(sim.is_in_bridge(sim.x_pos.at(0), sim.y_pos.at(0)));
    }


    BOOST_AUTO_TEST_CASE(test_explosion_bug2) {
        auto sim = get_sim(1);
        sim.gate_is_flat = true;
        sim.setup();
        sim.start(1);
        double delta = 0.001;
        const double pi = 3.141592653589793;
        sim.x_pos.at(0) = -sim.bridge_length / 2 - delta * delta;
        sim.y_pos.at(0) = sim.bridge_width / 2 - delta;
        sim.directions.at(0) = pi / 4;
        BOOST_CHECK(sim.is_going_in(0));
        BOOST_CHECK(not sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(not sim.is_in_bridge(sim.x_pos.at(0), sim.y_pos.at(0)));
        sim.compute_next_impact(0);
        sim.update(0);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(sim.is_in_bridge(sim.x_pos.at(0), sim.y_pos.at(0)));
        sim.update(0);
        BOOST_CHECK(sim.is_in_gate(sim.x_pos.at(0), sim.y_pos.at(0), sim.LEFT));
        BOOST_CHECK(sim.is_in_bridge(sim.x_pos.at(0), sim.y_pos.at(0)));
    }

    BOOST_AUTO_TEST_CASE(test_second_bridge) {
        // Test if the upper and lower bound of the bridge function works correctly
        // Todo: Try to identify the error that comes without the extension
    }

    BOOST_AUTO_TEST_CASE(test_periodic_boundary) {
        // Test from left to right, and from right to left, and see if angle and y are maintained as well
    }

    BOOST_AUTO_TEST_CASE(test_second_height_zero) {
        //Test if we correctly reduce to the first bridge
    }

    BOOST_AUTO_TEST_CASE(test_no_collision_for_urn_at_second_bridge) {
        // Test if in case of the second bridge no particles stick
    }


BOOST_AUTO_TEST_SUITE_END();