#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <unordered_map>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <random>
#include <cmath>

namespace Constants {
    // Fundamental constants
    extern const double Coulomb_Constant;
    extern const double Avogadro_Number;
    extern const double Boltzmann_Constant;
    extern const double Pi;

    // Activation Energies (in J)
    extern const double Ea_1;
    extern const double Ea_2;
    extern const double Ea_3;
    extern const double Ea_4;
    extern const double Ea_5;
    // Reaction enthalpy
    extern const double deltaH;

    // Simulation screen (in meters)
    extern const Eigen::Vector3d SCREEN;

    // Reaction coordinate (just a constant)
    extern const double dcoor;
    extern const double dt;
    // Atomic data
    std::unordered_map<std::string, int> get_atomic_radius_Covalent();
    std::unordered_map<std::string, int> get_atomic_radius_VanDerWaals();
    std::unordered_map<std::string, double> get_atomic_weight();
    std::mt19937& get_gen();
}

#endif // CONSTANTS_HPP