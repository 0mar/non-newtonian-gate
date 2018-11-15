//
// Created by Omar Richardson on 24/10/2018.
//

#include "simulation.h"

const double PI = 3.14159265358979323;
const double EPS = 1E-14;
#define px positions(particle,0)
#define py positions(particle,1)

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Simulation::Simulation(int num_particles, double gate_radius) : num_particles(num_particles) {
    // Make this const?
    // Make this static?
    rd = std::make_shared<std::random_device>();
    rng = std::make_shared<std::mt19937>((*rd)());
    unif_real = std::make_shared<std::uniform_real_distribution<double>>(0, 1);
    circle_radius = 1;
    circle_distance = 0.5;
    bridge_height = 0.3;
    this->gate_radius = gate_radius;
    gate_capacity = 3;
}

void Simulation::setup() {
    max_path = circle_distance + bridge_height + circle_radius * 4; // Upper bound for the longest path
    next_impact_times = Eigen::ArrayXd::Zero(num_particles);
    impact_times = Eigen::ArrayXd::Zero(num_particles);
    in_gate = Eigen::ArrayXi::Zero(num_particles);
    in_left = Eigen::ArrayXi::Ones(num_particles);
    in_right = Eigen::ArrayXi::Zero(num_particles);
    left_center_x = -circle_distance / 2 - circle_radius;
    right_center_x = circle_distance / 2 + circle_radius;
    positions = Eigen::ArrayXXd::Zero(num_particles, 2);
    next_positions = Eigen::ArrayXXd::Zero(num_particles, 2);
    directions = Eigen::ArrayXd::Zero(num_particles);
    next_directions = Eigen::ArrayXd::Zero(num_particles);
    couple_bridge();
}

void Simulation::start() {
    /**
     * Initiate all particles, to the left
     */
    time = 0;
    last_written_time = 0;
    double box_x_radius = circle_distance / 2 + circle_radius * 2;
    double box_y_radius = circle_radius;
    if (gate_radius >= box_x_radius) {
        throw std::invalid_argument("Gate radius too large; no initialization possible");
    }
    if (bridge_size / 2 > circle_radius) {
        throw std::invalid_argument("Bridge larger than circle; no initialization possible");
    }
    for (int particle = 0; particle < num_particles; particle++) {
        double pos_x = 0;
        double pos_y = 0;
        while (not is_in_left_circle(pos_x, pos_y) or is_in_gate_radius(pos_x, pos_y) or is_in_bridge(pos_x, pos_y)) {
            pos_x = ((*unif_real)(*rng) - 0.5) * box_x_radius * 2;
            pos_y = ((*unif_real)(*rng) - 0.5) * box_y_radius * 2;
        }
        positions(particle, 0) = pos_x;
        positions(particle, 1) = pos_y;
        directions(particle) = ((*unif_real)(*rng) - 0.5) * 2 * PI;
        in_left(particle) = 1;
        compute_next_impact(particle);
    }
    measure();
}

void Simulation::update(double write_dt) {
    // Todo: Suboptimal. If we use a heap, we don't need to loop.
    double next_impact = next_impact_times(0);
    int particle = 0;
    for (int p = 0; p < num_particles; p++) {
        if (next_impact > next_impact_times(p)) {
            next_impact = next_impact_times(p);
            particle = p;
        }
    }
    if (next_impact < time) {
        throw std::invalid_argument("Next impact in history, not possible");
    }
    std::cout << impact_times << std::endl;
    std::cout << next_impact_times << std::endl;
    std::cout << "Next impact at " << next_impact << " of particle " << particle << std::endl;
    while (next_impact > last_written_time) {
        write_positions_to_file(last_written_time + write_dt);
        last_written_time += write_dt;
        printf("Slicing on %.3f\n", last_written_time);
    }
    positions(particle, 0) = next_positions(particle, 0);
    positions(particle, 1) = next_positions(particle, 1);
    directions(particle) = next_directions(particle);
    impact_times(particle) = next_impact;
    time = next_impact;
    in_left(particle) = 0;
    in_right(particle) = 0;
    if (px < 0) {
        in_left(particle) = 1;
    } else if (px > 0) {
        in_right(particle) = 1;
    }
    if (is_in_gate_radius(px, py)) {
        if (in_gate(particle) == 0) {
            in_gate(particle) = 1;
//            printf("Num particles at time %.2f in gate is %d\n", time, in_gate.sum());
            currently_in_gate.push_back(particle);
            bool explodes = check_gate_explosion(); //  Todo: give offending particle as parameter
            if (explodes) {
                currently_in_gate.erase(std::remove(currently_in_gate.begin(), currently_in_gate.end(), particle),
                                        currently_in_gate.end());
                in_gate(particle) = 0;
            }
        }
    } else {
        if (in_gate(particle) == 1) {
            in_gate(particle) = 0;
            currently_in_gate.erase(std::remove(currently_in_gate.begin(), currently_in_gate.end(), particle),
                                    currently_in_gate.end());
        }
    }
    compute_next_impact(particle);
    measure();
}

bool Simulation::check_gate_explosion() {
    if (currently_in_gate.size() > gate_capacity) {
//        std::cout << "Explosion at " << time << std::endl;
        for (int particle: currently_in_gate) {
//            printf("Bouncing particle %d\n", particle);
            double x, y;
            get_current_position(particle, x, y);
//            printf("Position updated to (%.3f, %.3f)\n", x, y);
            if (not is_in_domain(x, y)) {
//                printf("Particle %d not in domain\n", particle);
            } else if (not is_in_gate_radius(x, y)) {
//                printf("Particle %d not in radius\n", particle);
//                printf("Distance from center: %.4f\n",sqrt(x*x + y*y));
            }
            px = x;
            py = y;
            directions(particle) = get_retraction_angle(particle);
            impact_times(particle) = time;
            compute_next_impact(particle);
//            printf("After boom, we get new positions at time %.2f\n",next_impact_times(particle));
        }
        return true;
    }
    return false;
}

void Simulation::measure() {
    measuring_times.push_back(time);
    total_left.push_back(in_left.sum());
    total_right.push_back(in_right.sum());
}

void Simulation::print_status() {
    printf("Time passed: %.2f\n", time);
    for (int particle = 0; particle < num_particles; particle++) {
        printf("Particle %d at \nPosition (%.4f, %.4f) at t=%.2f, angle %.2f pi\n", particle, px, py,
               impact_times(particle), directions(particle) / PI);
        printf("Planned impact at\nPosition (%.4f, %.4f) at t=%.2f, angle %.2f pi\n",
               next_positions(particle, 0),
               next_positions(particle, 1), next_impact_times(particle), next_directions(particle) / PI);
    }
    printf("Particles left: %d, particles right: %d\n", in_left.sum(), in_right.sum());
    printf("Particles in gate: %d\n", in_gate.sum());
}

void Simulation::write_positions_to_file(double time) {
    std::string filename = "results.dat";
    std::ofstream file;
    if (time == 0) {
        file.open(filename, std::ofstream::out | std::ofstream::trunc);
        file << "num_particles\tgate_radius\tcircle_radius\tcircle_distance\tbridge_height\tbridge_size\n";
        file << num_particles << " " << gate_radius << " " << circle_radius << " " << circle_distance << " "
             << bridge_height << " " << bridge_size << std::endl;
        file.close();
    }
    file.open(filename, std::ios_base::app);
    file << time << std::endl;
    for (int particle = 0; particle < num_particles; particle++) {
        file << px + (next_positions(particle, 0) - px) * (impact_times(particle) - time) /
                     (impact_times(particle) - next_impact_times(particle)) << " ";
    }
    file << std::endl;
    for (int particle = 0; particle < num_particles; particle++) {
        file << py + (next_positions(particle, 1) - py) * (impact_times(particle) - time) /
                     (impact_times(particle) - next_impact_times(particle)) << " ";
    }
    file << std::endl;
    for (int particle = 0; particle < num_particles; particle++) {
        file << directions(particle) << " ";
        printf("For particle %d, teller, %.3f, noemer: %.3f\n", particle, impact_times(particle) - time,
               impact_times(particle) - next_impact_times(particle));
    }
    file << std::endl;
    file.close();
}

void Simulation::write_totals_to_file() {
    std::string filename = "totals.dat";
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    for (double time: measuring_times) {
        file << time << "\t";
    }
    file << std::endl;
    for (int left: total_left) {
        file << left << "\t";
    }
    file << std::endl;
    for (int right: total_right) {
        file << right << "\t";
    }
    file << std::endl;
}

void Simulation::finish() {
    write_totals_to_file();
}

void Simulation::couple_bridge() {
    /**
     * A priori, the bridge does not connect to the circles.
     * We need to make the bridge a little bit longer so the ends connect too.
     * We do this by computing the intersections between the bridge lines and the circle
     */
    double lx = 0;
    double ly = bridge_height / 2;
    double angle;
    positions(0, 0) = lx;
    positions(0, 1) = ly;
    directions(0) = 0;
    bridge_size = 2 * this->time_to_hit_circle(0, right_center_x, angle);
}

bool Simulation::is_in_domain(double x, double y) {
    if (is_in_bridge(x, y)) {
        return true;
    } else {
        if (x < 0) {
            return is_in_left_circle(x, y);
        } else {
            return is_in_right_circle(x, y);
        }
    }
}

bool Simulation::is_in_left_circle(const double x, const double y) {
    return (x - left_center_x) * (x - left_center_x) + y * y < circle_radius * circle_radius;
}

bool Simulation::is_in_right_circle(const double x, const double y) {
    return (x - right_center_x) * (x - right_center_x) + y * y < circle_radius * circle_radius;
}

bool Simulation::is_in_bridge(const double x, const double y) {
/**
 * Note that these function is not mutually exclusive with left and right circle.
 */
    return abs(x) < bridge_size / 2 and abs(y) < bridge_height / 2;
}

bool Simulation::is_in_gate_radius(const double x, const double y) {
    /**
     * Check if the position is in the circular gate radius
     */
    return x * x + y * y < gate_radius * gate_radius;
}

void Simulation::compute_next_impact(const int particle) {
    /**
     * Start from some x, y, alpha.
     * Compute the location of boundary hit
     * Compute the time to next boundary hit
     * if boundary hit is gate hit:
     *  Flag this gate hit
     *
     *  How do we compute next boundary hit?
     *  We have 4 obstacles:
     *  Circle left, circle right, bridge, gate.
     *  Earliest impact count.
     *  If earliest impact is not a gate, then disregard this particle until it hits.
     *
     */
    double next_time = max_path;
    double next_angle = 0;
    double angle;
    double to_bridge = time_to_hit_bridge(particle, angle);
    // this flow is not supah dupah
    if (to_bridge < next_time) {
        next_time = to_bridge;
        next_angle = get_reflection_angle(directions(particle), angle);
    }
    double to_left = time_to_hit_circle(particle, left_center_x, angle);
    if (to_left < next_time) {
        next_time = to_left;
        next_angle = get_reflection_angle(directions(particle), angle);
    }
    double to_right = time_to_hit_circle(particle, right_center_x, angle);
    if (to_right < next_time) {
        next_time = to_right;
        next_angle = get_reflection_angle(directions(particle), angle);
    }
    double to_gate = time_to_hit_gate(particle);
    if (to_gate < next_time) {
        next_time = to_gate + EPS; // In the circle should be guaranteed in; out should be out
        next_angle = directions(particle);
        double nx = positions(particle, 0) + next_time * cos(directions(particle));
        double ny = positions(particle, 1) + next_time * sin(directions(particle));
        if (is_in_gate_radius(px, py) and is_in_gate_radius(nx, ny)) {
            printf("Small movement (%.3e) for particle %d detected\n", next_time, particle);
        }
    }
    if (angle == 0) {
        throw std::invalid_argument("wrong angle");
    }
    next_positions(particle, 0) = positions(particle, 0) + next_time * cos(directions(particle));
    next_positions(particle, 1) = positions(particle, 1) + next_time * sin(directions(particle));
    next_impact_times(particle) = time + next_time;
    next_directions(particle) = next_angle;
}

void Simulation::get_current_position(int particle, double &x, double &y) {
    /**
     * Interpolate position at the current time. Returns in referenced variables
     */
    x = px + (next_positions(particle, 0) - px) * (impact_times(particle) - time) /
             (impact_times(particle) - next_impact_times(particle));

    y = py + (next_positions(particle, 1) - py) * (impact_times(particle) - time) /
             (impact_times(particle) - next_impact_times(particle));
}

double Simulation::time_to_hit_bridge(const int particle, double &normal_angle) {
    /**
     * Check if we hit the bottom line, and check if we hit the top line, and return a float.
     */
    //Recall: px=positions(particle,0), py=positions(particle,1)
    double rx = max_path * cos(directions(particle));
    double ry = max_path * sin(directions(particle));
    double sx = bridge_size;
    double sy = 0;
    // q_bottom = (left_x, bottom_y) and q_top = (left_x, top_y)
    // u = (q − p) × r / (r × s)
    double u1 = ((-bridge_size / 2 - px) * ry - (-bridge_height / 2 - py) * rx) / (rx * sy - ry * sx);
    double u2 = ((-bridge_size / 2 - px) * ry - (bridge_height / 2 - py) * rx) / (rx * sy - ry * sx);
    // t = (q − p) × s / (r × s)
    double t1 = ((-bridge_size / 2 - px) * sy - (-bridge_height / 2 - py) * sx) / (rx * sy - ry * sx);
    double t2 = ((-bridge_size / 2 - px) * sy - (bridge_height / 2 - py) * sx) / (rx * sy - ry * sx); // can be faster
    double min_t = 1;
    if (EPS < t1 and t1 < min_t and 0 <= u1 and u1 <= 1) {
        min_t = t1 - EPS;
        normal_angle = PI / 2;
    }
    if (EPS < t2 and t2 < min_t and 0 <= u2 and u2 <= 1) {
        min_t = t2 - EPS;
        normal_angle = -PI / 2;
    }
    return min_t * max_path;
}

double Simulation::time_to_hit_circle(const int particle, const double center_x, double &normal_angle) {
    /**
     * Compute the time until next impact with one of the circle boundaries
     */
    double add_x = max_path * cos(directions(particle));
    double add_y = max_path * sin(directions(particle));
    const double t_pos_x = (px - center_x) / circle_radius;
    const double t_pos_y = (py - 0) / circle_radius;
    const double t_add_x = add_x / circle_radius;
    const double t_add_y = add_y / circle_radius;
    // Compose quadratic equation
    const double a = t_add_x * t_add_x + t_add_y * t_add_y;
    const double b = 2 * t_pos_x * t_add_x + 2 * t_pos_y * t_add_y;
    const double c = t_pos_x * t_pos_x + t_pos_y * t_pos_y - 1;
    const double D_sqrt = sqrt(b * b - 4 * a * c);
    if (D_sqrt < 0) printf("hmmmmm neg D\n");
    // compute roots
    const double t1 = (-b - D_sqrt) / (2 * a);
    const double t2 = (-b + D_sqrt) / (2 * a);
    // Find minimal root between 0 and 1 not in bridge
    double min_t = 1;
    double impact_x = 0;
    double impact_y = 0;
    if (EPS < t1 and t1 < min_t) {
        impact_x = px + t1 * add_x;
        impact_y = py + t1 * add_y;
        // Only hitting the circle if not in the bridge
        if (not is_in_bridge(impact_x, impact_y)) {
            normal_angle = atan2(0 - impact_y, center_x - impact_x);
            min_t = t1 - EPS;
        }
    }
    if (EPS < t2 and t2 < min_t) {
        impact_x = px + t2 * add_x;
        impact_y = py + t2 * add_y;
        // Only hitting the circle if not in the bridge
        if (not is_in_bridge(impact_x, impact_y)) {
            normal_angle = atan2(0 - impact_y, center_x - impact_x);
            min_t = t2 - EPS;
        }
    }
    return min_t * max_path;
}

double Simulation::get_reflection_angle(const double angle_in, const double normal_angle) {
    return fmod(2 * normal_angle - angle_in + PI, 2 * PI);
}

double Simulation::get_retraction_angle(const int particle) {
    int side = sgn(positions(particle, 0));
    double new_angle = ((*unif_real)(*rng) - 0.5) * PI + PI / 2 * (1 - sgn(side));
    return new_angle;
}

double Simulation::time_to_hit_gate(const int particle) {
    if (gate_radius > 0) {
        double add_x = max_path * cos(directions(particle));
        double add_y = max_path * sin(directions(particle));
        const double t_pos_x = px / gate_radius;
        const double t_pos_y = py / gate_radius;
        const double t_add_x = add_x / gate_radius;
        const double t_add_y = add_y / gate_radius;
        // Compose quadratic equation
        const double a = t_add_x * t_add_x + t_add_y * t_add_y;
        const double b = 2 * t_pos_x * t_add_x + 2 * t_pos_y * t_add_y;
        const double c = t_pos_x * t_pos_x + t_pos_y * t_pos_y - 1;
        const double D_sqrt = sqrt(b * b - 4 * a * c);
        if (D_sqrt < 0) printf("hmmmmm neg D\n");
        // compute roots
        const double t1 = (-b - D_sqrt) / (2 * a);
        const double t2 = (-b + D_sqrt) / (2 * a);
        // Find minimal root between 0 and 1
        double min_t = 1;
        if (EPS < t1 and t1 < min_t) {
            min_t = t1;
        }
        if (EPS < t2 and t2 < min_t) {
            min_t = t2;
        }
//    double impact_x = px + min_t * add_x;
//    double impact_y = py + min_t * add_y;
        return min_t * max_path;
    } else {
        return max_path;
    }

}
