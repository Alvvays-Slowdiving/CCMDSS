#ifndef FUNC_HPP
#define FUNC_HPP

#include "constants.hpp"
// 함수 선언
Eigen::Matrix3d Rotator_w(const Eigen::Vector3d& omega, double dt);
Eigen::Vector3d normalization(const Eigen::Vector3d& x);
Eigen::Vector3d E_field(const Eigen::Vector3d& dipole, const Eigen::Vector3d& r, double Coulomb_Constant);
double Potential_E(const Eigen::Vector3d& dipole, const Eigen::Vector3d& E);
Eigen::Vector3d Torque(const Eigen::Vector3d& dipole, const Eigen::Vector3d& E);
Eigen::Vector3d random_angular(std::mt19937& gen = Constants::get_gen());
double velocity_random(double m, double Temp, double Boltzmann_constant = Constants::Boltzmann_Constant, std::mt19937& gen = Constants::get_gen());
double uniform(std::mt19937& gen = Constants::get_gen(), double min = 0.0, double max = 1.0);
// double Q(const double& a, const double& x);
#endif // FUNC_H