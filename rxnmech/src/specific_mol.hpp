#ifndef SPECIFIC_MOL_HPP
#define SPECIFIC_MOL_HPP

#include "molecule.hpp"

typedef class NO2 : public molecule {
public:
    NO2();
    NO2(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius, const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM,
        const Eigen::Vector3d& torque, const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
}NO2;

typedef class CO : public molecule {
public:
    CO();
    CO(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius, const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM,
        const Eigen::Vector3d& torque, const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
}CO;

typedef class NO : public molecule {
public:
    NO();
    NO(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius, const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM,
        const Eigen::Vector3d& torque, const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
} NO;

typedef class NO3 : public molecule {
public:
    NO3();
    NO3(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius, const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM,
        const Eigen::Vector3d& torque, const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
} NO3;

typedef class Ar : public molecule {
public:
    Ar();
    Ar(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
} Ar;

typedef class Cl : public molecule {
public:
    Cl();
    Cl(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
} Cl;
typedef class CH4 : public molecule {
public:
    CH4();
    CH4(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius, \
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles = true);
} CH4;
vector<atom> NO2_initializer();
vector<atom> CO_initializer();
vector<atom> NO_initializer();
vector<atom> NO3_initializer();
vector<atom> Ar_initializer();
vector<atom> Cl_initializer();
vector<atom> CH4_initializer();
#endif // SPECIFIC_MOL_HPP