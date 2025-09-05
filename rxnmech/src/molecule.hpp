#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include "func.hpp"
#include "atom.hpp"
using namespace std;

typedef class molecule {
public:
    Eigen::Vector3d dipole;
    vector<atom> atoms;
    int numatoms;
    Eigen::Vector3d v_CM;
    Eigen::Vector3d coor_CM;
    Eigen::Vector3d w_CM;
    Eigen::Vector3d torque;
    Eigen::Vector3d force;
    vector<Eigen::Vector3d> atoms_coor;
    Eigen::Matrix3d inertia_tensor;
    double mass;
    bool is_polyatomic;
    bool is_linear;
    Eigen::Vector3d linear_axis;
    int degree_of_freedom;
    double radius;
    int index;
    string name;
    molecule(); // Default constructor
    molecule(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules,int num, double* max_radius, vector<atom> atoms_array, const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM\
        , const Eigen::Vector3d& torque, const Eigen::Vector3d& force, const Eigen::Vector3d& rotate, const Eigen::Vector3d& dipole = Eigen::Vector3d(0, 0, 0)\
        ,bool is_linear = false, bool save_particles = true); /*num은 원자 수인 것 같다.*/
    double Kinetic_Energy() const;
    double trans_Energy() const;
    double Rotation_Energy() const;
    void rotate(Eigen::Vector3d omega);
    void rotation(double dt);
    void update(vector<vector<vector<int>>>* checking_lst, int num_particles, double dt);
    bool check_collision(const molecule& other) const;
    bool collision(vector<vector<vector<int>>>* checking_lst, molecule& other, int* i, int* j);
    void cal_collision(molecule& other, const int& i, const int& j);
    bool check_collision_wall() const;
    void collision_wall();
    void cal_collision_wall(const int& i, const int& num);
    bool operator==(const molecule& rhs) const;
    bool operator!=(const molecule& rhs) const;
    void del(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules);
    Eigen::Vector3d get_w_CM(const double& Energy, const Eigen::Vector3d& vector) const;
    Eigen::Vector3d get_v_CM(const double& Energy, const Eigen::Vector3d& vector) const;
    virtual string get_name() const;
    virtual ~molecule(); // Default destructor
} molecule;
#endif // MOLECULE_HPP