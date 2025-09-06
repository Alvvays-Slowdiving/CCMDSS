#include "func.hpp"

using namespace Eigen;
using namespace std;

Matrix3d Rotator_w(const Vector3d& omega, double dt) {
    double norm_omega = omega.norm();
    Vector3d u = (norm_omega != 0) ? omega / norm_omega : Vector3d(1.0, 0.0, 0.0);
    double theta = norm_omega * dt;

    double ct = cos(theta);
    double st = sin(theta);
    double one_minus_ct = 1.0 - ct;

    Matrix3d R;
    R << ct + one_minus_ct * u(0) * u(0),        one_minus_ct * u(0) * u(1) - st * u(2), one_minus_ct * u(0) * u(2) + st * u(1),
         one_minus_ct * u(1) * u(0) + st * u(2), ct + one_minus_ct * u(1) * u(1),        one_minus_ct * u(1) * u(2) - st * u(0),
         one_minus_ct * u(2) * u(0) - st * u(1), one_minus_ct * u(2) * u(1) + st * u(0), ct + one_minus_ct * u(2) * u(2);

    return R;
}
Vector3d normalization(const Vector3d& x) {
    double norm_x = x.norm();
    return (norm_x != 0) ? x / norm_x : Vector3d(1.0, 0.0, 0.0);
}
Vector3d E_field(const Vector3d& dipole, const Vector3d& r, double Coulomb_Constant) {
    Vector3d r_hat = r.normalized();
    double r_norm = r.norm();
    return (Coulomb_Constant / pow(r_norm, 3)) * (3.0 * dipole.dot(r_hat) * r_hat - dipole);
}
double Potential_E(const Vector3d& dipole, const Vector3d& E) {
    return -dipole.dot(E);
}
Vector3d Torque(const Vector3d& dipole, const Vector3d& E) {
    return dipole.cross(E);
}
Vector3d random_angular(std::mt19937& gen) {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    std::uniform_real_distribution<> angle_dist(0.0, 2 * 3.14159265358979323846); // 2 * pi

    Vector3d v(dist(gen), dist(gen), dist(gen));
    return normalization(v) * angle_dist(gen);
}
double velocity_random(double m, double Temp, double Boltzmann_constant, std::mt19937& gen) {
    double scale = sqrt(Boltzmann_constant * Temp / m);
    std::normal_distribution<> normal(0.0, scale);
    // For a 3D speed, you usually sample 3 normal distributions and take norm
    Vector3d v(normal(gen), normal(gen), normal(gen));
    return v.norm();  // Equivalent to maxwell.rvs() in 3D
}
double uniform(std::mt19937& gen, double min, double max) {
    std::uniform_real_distribution<> dist(min, max);
    return dist(gen);
}
// double Q(const double& a, const double& x) {
//     // Q function implementation using GSL
//     return gsl_sf_gamma_inc_Q(a, x);
// }