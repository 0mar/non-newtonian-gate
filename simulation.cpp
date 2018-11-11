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
    bridge_height = 0.1;
    this->gate_radius = gate_radius;
    gate_capacity = 1;
}

void Simulation::setup() {
    max_path = circle_distance + bridge_height + circle_radius * 4; // Upper bound for the longest path
    next_impact_times = Eigen::ArrayXd::Zero(num_particles); // Todo: This should be a heap
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

void Simulation::update() {
    // Todo: Suboptimal. If we use a heap, we don't need to loop.
    double next_impact = next_impact_times(0);
    int particle = 0;
    for (int p = 0; p < num_particles; p++) {
        if (next_impact > next_impact_times(p)) {
            next_impact = next_impact_times(p);
            particle = p;
        }
    }
    impact_times(particle) = time;
    positions(particle, 0) = next_positions(particle, 0);
    positions(particle, 1) = next_positions(particle, 1);
    directions(particle) = next_directions(particle);
    in_left(particle) = 0;
    in_right(particle) = 0;
    if (is_in_left_circle(px, py)) {
        in_left(particle) = 1;
    } else if (is_in_right_circle(px, py)) {
        in_right(particle) = 1;
    }
    if (is_in_gate_radius(px, py)) {
        in_gate(particle) = 1;
        check_gate_explosion();
    } else {
        in_gate(particle) = 0;
    }
    time = next_impact;
    compute_next_impact(particle);
    measure();
}

void Simulation::check_gate_explosion() {
    if (in_gate.sum() > gate_capacity) {
        // printf("Explosion! Over %d particles in the gate\n", gate_capacity);
    }
}

void Simulation::measure() {
    measuring_times.push_back(time);
    total_left.push_back(in_left.sum());
    total_right.push_back(in_right.sum());
}

void Simulation::print_status() {
    printf("Time passed: %.2f\n", time);
    for (int particle = 0; particle < num_particles; particle++) {
        printf("Particle %d at position (%.2f, %.2f) at Dt %.2f, angle %.2f pi\n", particle, px, py,
               next_impact_times(particle), directions(particle) / PI);
    }
    printf("Particles left: %d, particles right: %d\n", in_left.sum(), in_right.sum());
    printf("Particles in gate: %d\n", in_gate.sum());
}

void Simulation::write_to_file(bool interpolate) {
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
    if (interpolate) {
        for (int particle = 0; particle < num_particles; particle++) {
            file << px + (next_positions(particle, 0) - px) * (impact_times(particle) - time) /
                         (impact_times(particle) - next_impact_times(particle)) << " ";
        }
        file << std::endl;
        for (int particle = 0; particle < num_particles; particle++) {
            file << py + (next_positions(particle, 1) - py) * (impact_times(particle) - time) /
                         (impact_times(particle) - next_impact_times(particle)) << " ";
        }
    } else {
        for (int particle = 0; particle < num_particles; particle++) {
            file << px << " ";
        }
        file << std::endl;
        for (int particle = 0; particle < num_particles; particle++) {
            file << py << " ";
        }
    }
    file << std::endl;
    for (int particle = 0; particle < num_particles; particle++) {
        file << directions(particle) << " ";
    }
    file << std::endl;
    file.close();
}

void Simulation::finish() {
    for (unsigned long i = 0; i < measuring_times.size(); i++) {
        printf("Time: %.2f\tLeft:%d\tRight:%d\n", measuring_times.at(i), total_left.at(i), total_right.at(i));
    }
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
    if (to_gate < next_time) { // keep last: we change stuff here.
        next_time = to_gate;
        next_angle = directions(particle);
    }
    next_positions(particle, 0) = positions(particle, 0) + next_time * cos(directions(particle));
    next_positions(particle, 1) = positions(particle, 1) + next_time * sin(directions(particle));
    next_impact_times(particle) = time + next_time;
    next_directions(particle) = next_angle;
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
        // Find minimal root between 0 and 1 not in bridge
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
