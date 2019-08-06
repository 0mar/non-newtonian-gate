//
// Created by Omar Richardson on 24/10/2018.
// Hints for optimization:
// 1. Be careful with compute_next_impact() calls
//

#include "simulation.h"

const double PI = 3.14159265358979323;
const double EPS = 1E-14;
#define px x_pos[particle]
#define py y_pos[particle]

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Simulation::Simulation(const int &num_particles, const double &gate_radius) : num_particles(num_particles),
                                                                              gate_radius(gate_radius) {
    rd = std::make_shared<std::random_device>();
    rng = std::make_shared<std::mt19937>((*rd)());
    unif_real = std::make_shared<std::uniform_real_distribution<double>>(0, 1);
    circle_radius = 1;
    circle_distance = 0.5;
    bridge_height = 0.3;
    bridge_size = 0; // Computed later
    left_gate_capacity = 3;
    right_gate_capacity = 3;
    explosion_direction_is_random = true;
}

void Simulation::setup() {
    max_path = circle_distance + bridge_height + circle_radius * 4; // Upper bound for the longest path
    next_impact_times.reserve(num_particles); // Will be filled later
    next_impact_times_check.resize(num_particles); // Will be filled later
    last_written_times.resize(num_particles);
    impact_times.resize(num_particles);
    in_left_gate.resize(num_particles);
    in_right_gate.resize(num_particles);
    x_pos.resize(num_particles);
    y_pos.resize(num_particles);
    next_x_pos.resize(num_particles);
    next_y_pos.resize(num_particles);
    directions.resize(num_particles);
    next_directions.resize(num_particles);
    gate_arrays.push_back(in_left_gate);
    gate_arrays.push_back(in_right_gate);
    gate_contents.push_back(currently_in_left_gate);
    gate_contents.push_back(currently_in_right_gate);
    gate_capacities.push_back(left_gate_capacity);
    gate_capacities.push_back(right_gate_capacity);

    left_center_x = -circle_distance / 2 - circle_radius;
    right_center_x = circle_distance / 2 + circle_radius;

    couple_bridge();
}

void Simulation::start() {
    /**
     * Initiate all particles, to the left
     */
    time = 0;
    last_written_time = 0;
    in_left = 0;
    in_right = 0;
    const double box_x_radius = circle_distance / 2 + circle_radius * 2;
    const double box_y_radius = circle_radius;
    if (gate_radius >= box_x_radius) {
        throw std::invalid_argument("Gate radius too large; no initialization possible");
    }
    if (bridge_size / 2 > circle_radius) {
        throw std::invalid_argument("Bridge larger than circle; no initialization possible");
    }
    for (unsigned long particle = 0; particle < num_particles; particle++) {
        double pos_x = 0;
        double pos_y = 0;
        while (not is_in_circle(pos_x, pos_y, LEFT) or is_in_gate(pos_x, pos_y, LEFT) or
               is_in_bridge(pos_x, pos_y)) {
            pos_x = ((*unif_real)(*rng) - 0.5) * box_x_radius * 2;
            pos_y = ((*unif_real)(*rng) - 0.5) * box_y_radius * 2;
        }
        px = pos_x;
        py = pos_y;
        directions.at(particle) = ((*unif_real)(*rng) - 0.5) * 2 * PI;
        compute_next_impact(particle);
        in_left++;
    }
    measure();
}

void Simulation::start_evenly() {
    /**
     * Initiate all particles, as many left as right
     */
    time = 0;
    last_written_time = 0;
    in_left = 0;
    in_right = 0;
    const double box_x_radius = circle_distance / 2 + circle_radius * 2;
    const double box_y_radius = circle_radius;
    if (gate_radius >= box_x_radius) {
        throw std::invalid_argument("Gate radius too large; no initialization possible");
    }
    if (bridge_size / 2 > circle_radius) {
        throw std::invalid_argument("Bridge larger than circle; no initialization possible");
    }
    for (unsigned long particle = 0; particle < num_particles / 2; particle++) {
        double pos_x = 0;
        double pos_y = 0;
        while (not is_in_circle(pos_x, pos_y, LEFT) or is_in_gate(pos_x, pos_y, LEFT) or
               is_in_bridge(pos_x, pos_y)) {
            pos_x = ((*unif_real)(*rng) - 0.5) * box_x_radius * 2;
            pos_y = ((*unif_real)(*rng) - 0.5) * box_y_radius * 2;
        }
        px = pos_x;
        py = pos_y;
        directions.at(particle) = ((*unif_real)(*rng) - 0.5) * 2 * PI;
        compute_next_impact(particle);
        in_left++;
    }
    for (unsigned long particle = (unsigned long) num_particles / 2; particle < num_particles; particle++) {
        double pos_x = 0;
        double pos_y = 0;
        while (not is_in_circle(pos_x, pos_y, RIGHT) or is_in_gate(pos_x, pos_y, RIGHT) or
               is_in_bridge(pos_x, pos_y)) {
            pos_x = ((*unif_real)(*rng) - 0.5) * box_x_radius * 2;
            pos_y = ((*unif_real)(*rng) - 0.5) * box_y_radius * 2;
        }
        px = pos_x;
        py = pos_y;
        directions.at(particle) = ((*unif_real)(*rng) - 0.5) * 2 * PI;
        compute_next_impact(particle);
        in_right++;
    }
    measure();
}

void Simulation::get_next_impact(unsigned long &particle, double &next_impact) {
    bool impact_found = false;
    double next_impact_check = next_impact_times_check.at(0);

    for (const Impact &impact: next_impact_times) {
        printf("Particle number: %lu at time %.3f (%.3f) - %.2f\n",impact.particle,impact.time,next_impact_times_check[impact.particle],impact.stamp);
    }
    unsigned long particle_check = 0;
    for (unsigned long p = 0; p < num_particles; p++) {
        if (next_impact_check > next_impact_times_check.at(p)) {
            next_impact_check = next_impact_times_check.at(p);
            particle_check = p;
        }
    }

    while (not impact_found and not next_impact_times.empty()) {
        Impact &impact = next_impact_times.front();
//        printf("Candidate impact: particle: %lu with time %.2f\n",impact.particle,impact.time);
        if (impact.stamp == last_written_times[impact.particle]) {
            // This only works because we are not doing any floating point arithmetic
            impact_found = true;
            particle = impact.particle;
            next_impact = impact.time;
        } else {
            printf("Skipping ")
            std::cout << "Skipping " << impact.particle << std::endl;
        }
        std::pop_heap(next_impact_times.begin(), next_impact_times.end());
        next_impact_times.pop_back();
    }
    if (not impact_found) {
        throw std::invalid_argument("No impact found with a correct time stamp");
    }

    if  (not (particle==particle_check and next_impact == next_impact_check)) {
        printf("\nError in computing min distance: particle %lu (%lu), time diff %.2e \n",particle,particle_check,next_impact_check - next_impact);
        printf("%.2f\n",time);
        exit(1);
    }
    if (time>0.01) {

    }
//    printf("Decided impact: particle %lu with time %.2f\n",particle,next_impact);
}

void Simulation::update(const double &write_dt) {
    // Find next event: the first particle that has a new impact
    unsigned long particle;
    double next_impact;
    get_next_impact(particle, next_impact);
//    printf("Event list size: %lu/%d\n",next_impact_times.size(),num_particles);
    // Write a time slice, if desired
    if (write_dt > 0) {
        while (next_impact > last_written_time + write_dt) {
            write_positions_to_file(last_written_time + write_dt);
            last_written_time += write_dt;
        }
    }

    // Update the data of the particle with the collision
    if (not is_in_domain(next_x_pos[particle], next_y_pos[particle])) {
        // printf("Stray particle %d about to leave domain at (%.5f,%.5f), re-entered\n", (int) particle, px, py);
        next_x_pos[particle] = sgn(next_x_pos) * (circle_distance / 2 + circle_radius);
        next_y_pos[particle] = 0;
    }

    // Process the location of the particle
    if (px > 0 and next_x_pos[particle] < 0) {
        in_right--;
        in_left++;
    } else if (px < 0 and next_x_pos[particle] > 0) {
        in_left--;
        in_right++;
    } else if (px == 0) {
        std::cout << "Exactly zero position, administration just lost a particle" << std::endl;
    }
    px = next_x_pos[particle];
    py = next_y_pos[particle];
    directions[particle] = next_directions[particle];
    impact_times[particle] = next_impact;
    time = next_impact;


    // Check if the particle explodes
    for (unsigned long direction = 0; direction < 2; direction++) {
        if (is_in_gate(px, py, direction)) {
            check_gate_admission(particle, direction);
        } else {
            check_gate_departure(particle, direction);
        }
    }

    // Find out when the next collision takes place
    compute_next_impact(particle);

    // Do something useful with this information
    measure();
}

bool Simulation::is_in_gate(const double &x, const double &y, const unsigned long &direction) const {
    return ((int) direction * 2 - 1) * x >= 0 and x * x + y * y < gate_radius * gate_radius;
}

void Simulation::check_gate_admission(const unsigned long &particle, const unsigned long &direction) {
    if (gate_arrays[direction][particle] == 0) {
        // Not yet in gate, check admission
        if (gate_contents[direction].size() > gate_capacities[direction] - 1) {
            explode_gate(particle, direction);
        } else {
            gate_contents[direction].push_back(particle);
            gate_arrays[direction][particle] = 1;
        }
    }
}

void Simulation::check_gate_departure(const unsigned long &particle, const unsigned long &direction) {
    if (gate_arrays[direction][particle] == 1) {
        // Freshly leaving the gate
        gate_contents[direction].erase(std::remove(gate_contents[direction].begin(),
                                                   gate_contents[direction].end(), particle),
                                       gate_contents[direction].end());
        gate_arrays[direction][particle] = 0;
    }
}

void Simulation::explode_gate(const unsigned long &exp_particle, const unsigned long &direction) {
    do {
        directions[exp_particle] = get_retraction_angle(exp_particle);
        compute_next_impact(exp_particle);
    } while (not is_in_domain(next_x_pos[exp_particle], next_y_pos[exp_particle]));
    for (unsigned long particle: gate_contents[direction]) {
//            printf("Bouncing particle %d with position (%.3f, %.3f)\n", particle,px,py);
//            printf("Last impact time: %.2f\tNext impact time%.2f\n",impact_times(particle),next_impact_times(particle));
        double x, y;
        get_current_position(particle, x, y);
//            printf("Position updated to (%.3f, %.3f)\n", x, y);
        if (not is_in_domain(x, y)) {
            printf("Particle %d not in domain\n", (int) particle);
        } else if (not is_in_gate(x, y, direction)) {
            // printf("Particle %d not in radius, distance from center: %.4f. Removed from gate list\n", (int) particle, sqrt(x * x + y * y));
            gate_contents[direction].erase(std::remove(gate_contents[direction].begin(),
                                                       gate_contents[direction].end(), particle),
                                           gate_contents[direction].end()); //Todo: Is removing during loop smart?
        }
        px = x;
        py = y;
        directions[particle] = get_retraction_angle(particle);
        impact_times[particle] = time;
        compute_next_impact(particle);
//            printf("After boom, we get new positions at time %.2f\n",next_impact_times(particle));
    }
}

void Simulation::measure() {
    measuring_times.push_back(time);
    total_left.push_back(in_left);
}

void Simulation::print_status() const {
    printf("Time passed: %.2f\n", time);
    for (unsigned long particle = 0; particle < num_particles; particle++) {
        printf("Particle %d at \nPosition (%.4f, %.4f) at t=%.2f, angle %.2f pi\n", (int) particle, px, py,
               impact_times[particle], directions[particle] / PI);
        printf("Planned impact at\nPosition (%.4f, %.4f) at t=%.2f, angle %.2f pi\n",
               next_x_pos[particle],
               next_y_pos[particle], next_impact_times[particle].time, next_directions[particle] / PI);
    }
    printf("Particles left: %d, particles right: %d\n", (int) in_left, (int) in_right);
    printf("Particles in left gate: %d\t in right gate %d\n", (int) currently_in_left_gate.size(),
           (int) currently_in_right_gate.size());
}

void Simulation::write_positions_to_file(const double &time) const {
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
    for (unsigned long particle = 0; particle < num_particles; particle++) {
        file << px + (next_x_pos[particle] - px) * (impact_times[particle] - time) /
                     (impact_times[particle] - next_impact_times[particle].time) << " ";
    }
    file << std::endl;
    for (unsigned long particle = 0; particle < num_particles; particle++) {
        file << py + (next_y_pos[particle] - py) * (impact_times[particle] - time) /
                     (impact_times[particle] - next_impact_times[particle].time) << " ";
    }
    file << std::endl;
    for (unsigned long particle = 0; particle < num_particles; particle++) {
        file << directions[particle] << " ";
    }
    file << std::endl;
    file.close();
}

void Simulation::write_totals_to_file() const {
    std::string filename = "totals.dat";
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    for (double m_time: measuring_times) {
        file << m_time << "\t";
    }
    file << std::endl;
    for (unsigned long left: total_left) {
        file << left << "\t";
    }
    file << std::endl;
    for (unsigned long left: total_left) {
        file << num_particles - left << "\t";
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
    x_pos.at(0) = lx;
    y_pos.at(0) = ly;
    directions.at(0) = 0;
    bridge_size = 2 * this->time_to_hit_circle(0, right_center_x, angle);
}

bool Simulation::is_in_domain(const double &x, const double &y) const {
    if (is_in_bridge(x, y)) {
        return true;
    } else {
        if (x < 0) {
            return is_in_circle(x, y, LEFT);
        } else {
            return is_in_circle(x, y, RIGHT);
        }
    }
}

bool Simulation::is_in_circle(const double &x, const double &y, const unsigned long &side) const {
    if (side == LEFT) {
        return (x - left_center_x) * (x - left_center_x) + y * y < circle_radius * circle_radius;
    } else {
        return (x - right_center_x) * (x - right_center_x) + y * y < circle_radius * circle_radius;
    }
}


bool Simulation::is_in_bridge(const double &x, const double &y) const {
/**
 * Note that these function is not mutually exclusive with left and right circle.
 */
    return std::abs(x) < bridge_size / 2 and std::abs(y) < bridge_height / 2;
}

void Simulation::compute_next_impact(const unsigned long &particle) {
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
        next_angle = get_reflection_angle(directions[particle], angle);
    }
    double to_left = time_to_hit_circle(particle, left_center_x, angle);
    if (to_left < next_time) {
        next_time = to_left;
        next_angle = get_reflection_angle(directions[particle], angle);
    }
    double to_right = time_to_hit_circle(particle, right_center_x, angle);
    if (to_right < next_time) {
        next_time = to_right;
        next_angle = get_reflection_angle(directions[particle], angle);
    }
    double to_gate = time_to_hit_gate(particle);
    if (to_gate < next_time) {
        next_time = to_gate + EPS; // In the circle should be guaranteed in; out should be out
        next_angle = directions[particle];
//        if (is_in_gate_radius(px, py) and is_in_gate_radius(nx, ny)) {
//            printf("Small movement (%.3e) for particle %d detected\n", next_time, particle);
//        }
    }
    double to_middle = time_to_hit_middle(particle);

    if (to_middle < next_time) {
        next_time = to_middle + EPS;
        next_angle = directions[particle];
    }
    if (next_time == max_path) {
        printf("Next time, %.2f, maxpath %.2f\n", next_time, max_path);
        throw std::invalid_argument("No collision found");

    }
    next_x_pos[particle] = px + next_time * cos(directions[particle]);
    next_y_pos[particle] = py + next_time * sin(directions[particle]);
    unsigned long p = particle;
    next_impact_times_check[p] = time+next_time;
    next_impact_times.emplace_back(p, time + next_time, time);
    last_written_times[p] = time;
    std::push_heap(next_impact_times.begin(), next_impact_times.end());
    next_directions[particle] = next_angle;
}

void Simulation::get_current_position(const unsigned long &particle, double &x, double &y) const {
    /**
     * Interpolate position at the current time. Returns in referenced variables
     */
    if (impact_times[particle] == next_impact_times[particle].time) {
        x = px;
        y = py;
    } else {
        x = px + (next_x_pos[particle] - px) * (impact_times[particle] - time) /
                 (impact_times[particle] - next_impact_times[particle].time);
        y = py + (next_y_pos[particle] - py) * (impact_times[particle] - time) /
                 (impact_times[particle] - next_impact_times[particle].time);
    }
}

double Simulation::time_to_hit_bridge(const unsigned long &particle, double &normal_angle) const {
    /**
     * Check if we hit the bottom line, and check if we hit the top line, and return a float.
     */
    //Recall: px=positions(particle,0), py=positions(particle,1)
    double rx = max_path * cos(directions[particle]);
    double ry = max_path * sin(directions[particle]);
    double sx = bridge_size;
    double sy = 0;
    // q_bottom = (left_x, bottom_y) and q_top = (left_x, top_y)
    // u = (q − p) × r / (r × s)
    double u1 = ((-bridge_size / 2 - px) * ry - (-bridge_height / 2 - py) * rx) / (rx * sy - ry * sx);
    double u2 = ((-bridge_size / 2 - px) * ry - (bridge_height / 2 - py) * rx) / (rx * sy - ry * sx);
    // t = (q − p) × s / (r × s)
    double t1 = ((-bridge_size / 2 - px) * sy - (-bridge_height / 2 - py) * sx) / (rx * sy - ry * sx);
    double t2 = ((-bridge_size / 2 - px) * sy - (bridge_height / 2 - py) * sx) / (rx * sy - ry * sx); //could be faster
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

double
Simulation::time_to_hit_circle(const unsigned long &particle, const double &center_x, double &normal_angle) const {
    /**
     * Compute the time until next impact with one of the circle boundaries
     */
    double min_t = 1;
    double add_x = max_path * cos(directions[particle]);
    double add_y = max_path * sin(directions[particle]);
    const double t_pos_x = (px - center_x) / circle_radius;
    const double t_pos_y = (py - 0) / circle_radius;
    const double t_add_x = add_x / circle_radius;
    const double t_add_y = add_y / circle_radius;
    // Compose quadratic equation
    const double a = t_add_x * t_add_x + t_add_y * t_add_y;
    const double b = 2 * t_pos_x * t_add_x + 2 * t_pos_y * t_add_y;
    const double c = t_pos_x * t_pos_x + t_pos_y * t_pos_y - 1;
    const double D = b * b - 4 * a * c;
    if (D < 0) {
        min_t = 1;
    } else {
        // compute roots
        const double t1 = (-b - sqrt(D)) / (2 * a);
        const double t2 = (-b + sqrt(D)) / (2 * a);
        // Find minimal root between 0 and 1 not in bridge
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
    }
    return min_t * max_path;
}

double Simulation::get_reflection_angle(const double &angle_in, const double &normal_angle) const {
    return fmod(2 * normal_angle - angle_in + PI, 2 * PI);
}

double Simulation::get_retraction_angle(const unsigned long &particle) const {
    if (explosion_direction_is_random) {
        int side = sgn(px);
        return ((*unif_real)(*rng) - 0.5) * PI + PI / 2 * (1 - sgn(side));
    } else {
        if (cos(directions[particle]) * x_pos[particle] < 0) {
            return directions[particle] + PI;
        } else {
            return directions[particle];
        }
    }
}

double Simulation::time_to_hit_gate(const unsigned long &particle) const {
    /**
     * Compute time towards the circular gate by transforming the domain and solving a quadratic equation.
     * In order to ensure positivity of the solution, we use numerical rounding with epsilon for finding roots
     */
    if (gate_radius > 0) {
        double add_x = max_path * cos(directions[particle]);
        double add_y = max_path * sin(directions[particle]);
        const double t_pos_x = px / gate_radius;
        const double t_pos_y = py / gate_radius;
        const double t_add_x = add_x / gate_radius;
        const double t_add_y = add_y / gate_radius;
        // Compose quadratic equation
        const double a = t_add_x * t_add_x + t_add_y * t_add_y;
        const double b = 2 * t_pos_x * t_add_x + 2 * t_pos_y * t_add_y;
        const double c = t_pos_x * t_pos_x + t_pos_y * t_pos_y - 1;
        // Assume that D is indeed positive because we do not want to pay the price and check
        const double D_sqrt = sqrt(b * b - 4 * a * c);
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

double Simulation::time_to_hit_middle(const unsigned long &particle) const {
    /**
     * Uses a line-line intersection algorithm (with identical nomenclature) from 
     * https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
     */
    double min_t = 1;
    if (gate_radius > 0) {
        double rx = max_path * cos(directions[particle]);
        double ry = max_path * sin(directions[particle]);
        double sx = 0;
        double sy = gate_radius * 2;
        // u = (q − p) × r / (r × s)
        double denum = 1. / (rx * sy - ry * sx);
        double u = ((0 - px) * ry - (-gate_radius - py) * rx) * denum;
        // t = (q − p) × s / (r × s)
        double t = ((0 - px) * sy - (-gate_radius - py) * sx) * denum;
        if (EPS < t and t < min_t and 0 <= u and u <= 1) {
            min_t = t + EPS;
        }
    }
    return min_t * max_path;
}
