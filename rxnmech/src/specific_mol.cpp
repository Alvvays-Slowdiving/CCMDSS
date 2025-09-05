#include "specific_mol.hpp"

vector<atom> NO2_initializer() {
    vector<atom> NO2 = vector<atom>(3);
    NO2[0] = atom(Constants::get_atomic_radius_Covalent().at("N") * 1.0e-12, 
                   Constants::get_atomic_weight().at("N") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "N");
    NO2[1] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(110.3, 0.0, 46.482) * 1.0e-12, "O");
    NO2[2] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(-110.3, 0.0, 46.482) * 1.0e-12, "O");
    return NO2;
}
vector<atom> CO_initializer() {
    vector<atom> CO = vector<atom>(2);
    CO[0] = atom(Constants::get_atomic_radius_Covalent().at("C") * 1.0e-12, 
                   Constants::get_atomic_weight().at("C") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "C");
    CO[1] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 112.8) * 1.0e-12, "O");
    return CO;
}
vector<atom> NO_initializer() {
    vector<atom> NO = vector<atom>(2);
    NO[0] = atom(Constants::get_atomic_radius_Covalent().at("N") * 1.0e-12, 
                   Constants::get_atomic_weight().at("N") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "N");
    NO[1] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 115.0) * 1.0e-12, "O");
    return NO;
}
vector<atom> NO3_initializer() {
    vector<atom> NO3 = vector<atom>(4);
    NO3[0] = atom(Constants::get_atomic_radius_Covalent().at("N") * 1.0e-12, 
                   Constants::get_atomic_weight().at("N") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "N");
    NO3[1] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 123.9, 0.0) * 1.0e-12, "O");
    NO3[2] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(sqrt(3)*123.9/2, -123.9/2, 0.0) * 1.0e-12, "O");
    NO3[3] = atom(Constants::get_atomic_radius_Covalent().at("O") * 1.0e-12, 
                   Constants::get_atomic_weight().at("O") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(-sqrt(3)*123.9/2, -123.9/2, 0.0) * 1.0e-12, "O");
    return NO3;
}
vector<atom> Ar_initializer() {
    vector<atom> Ar = vector<atom>(1);
    Ar[0] = atom(Constants::get_atomic_radius_Covalent().at("Ar") * 1.0e-12, 
                 Constants::get_atomic_weight().at("Ar") / (Constants::Avogadro_Number * 1e+3), 
                 Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "Ar");
    return Ar;
}
vector<atom> Cl_initializer() {
    vector<atom> Cl = vector<atom>(1);
    Cl[0] = atom(Constants::get_atomic_radius_Covalent().at("Cl") * 1.0e-12, 
                 Constants::get_atomic_weight().at("Cl") / (Constants::Avogadro_Number * 1e+3), 
                 Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "Cl");
    return Cl;
}
vector<atom> CH4_initializer() {
    vector<atom> CH4 = vector<atom>(5);
    CH4[0] = atom(Constants::get_atomic_radius_Covalent().at("C") * 1.0e-12, 
                   Constants::get_atomic_weight().at("C") / (Constants::Avogadro_Number * 1e+3), 
                   Eigen::Vector3d(0.0, 0.0, 0.0) * 1.0e-12, "C");
    CH4[1] = atom(Constants::get_atomic_radius_Covalent().at("H") * 1.0e-12,
                   Constants::get_atomic_weight().at("H") / (Constants::Avogadro_Number * 1e+3),
                   Eigen::Vector3d(0.0, 0.0, 109.1) * 1.0e-12, "H");
    CH4[2] = atom(Constants::get_atomic_radius_Covalent().at("H") * 1.0e-12,
                    Constants::get_atomic_weight().at("H") / (Constants::Avogadro_Number * 1e+3),
                    Eigen::Vector3d(89.0, -51.4, -36.4) * 1.0e-12, "H"); 
    CH4[3] = atom(Constants::get_atomic_radius_Covalent().at("H") * 1.0e-12,
                    Constants::get_atomic_weight().at("H") / (Constants::Avogadro_Number * 1e+3),
                    Eigen::Vector3d(-89.0, -51.4, -36.4) * 1.0e-12, "H");
    CH4[4] = atom(Constants::get_atomic_radius_Covalent().at("H") * 1.0e-12,
                    Constants::get_atomic_weight().at("H") / (Constants::Avogadro_Number * 1e+3),
                    Eigen::Vector3d(0.0, 102.8, -36.4) * 1.0e-12, "H");
    return CH4;
}
NO2::NO2() : molecule() {
    this->numatoms = 3;
    this->atoms = NO2_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, 1.067405105e-30);
    this->is_polyatomic = true;
    this->is_linear = false;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "NO2";
}
NO2::NO2(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 3, max_radius, NO2_initializer(), v_CM, w_CM, torque, force, rotate, dipole = Eigen::Vector3d(0.0,0.0,1.067405105e-30), is_linear = false, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
        this->name = "NO2";
    }
CO::CO() : molecule() {
    this->numatoms = 2;
    this->atoms = CO_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, -4.069481961e-31);
    this->is_polyatomic = true;
    this->is_linear = true;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "CO";
}
CO::CO(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 2, max_radius, CO_initializer(), v_CM, w_CM, torque, force, rotate, dipole = Eigen::Vector3d(0.0,0.0,-4.069481961e-31), is_linear = true, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
        this->name = "CO";
    }
NO::NO() : molecule() {
    this->numatoms = 2;
    this->atoms = NO_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, -5.303669114e-31);
    this->is_polyatomic = true;
    this->is_linear = true;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "NO";
}
NO::NO(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 2, max_radius, NO_initializer(), v_CM, w_CM, torque, force, rotate, dipole = Eigen::Vector3d(0.0,0.0,-5.303669114e-31), is_linear = true, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
        this->name = "NO";
    }
NO3::NO3() : molecule() {
    this->numatoms = 4;
    this->atoms = NO3_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, 1.067405105e-30);
    this->is_polyatomic = true;
    this->is_linear = false;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "NO3";
}
NO3::NO3(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 4, max_radius, NO3_initializer(), v_CM, w_CM, torque, force, rotate, dipole = Eigen::Vector3d(0.0,0.0,0.0), is_linear = false, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
        this->name = "NO3";
    }
Ar::Ar() : molecule() {
    this->numatoms = 1;
    this->atoms = Ar_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->is_polyatomic = false;
    this->is_linear = false;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "Ar";
}
Ar::Ar(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 1, max_radius, Ar_initializer(), v_CM, Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0), force, Eigen::Vector3d(0.0, 0.0, 0.0), dipole = Eigen::Vector3d(0.0, 0.0, 0.0), is_linear = false, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        this->atoms[0].coor = this->coor_CM;
        this->name = "Ar";
    }
Cl::Cl() : molecule() {
    this->numatoms = 1;
    this->atoms = Cl_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->is_polyatomic = false;
    this->is_linear = false;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "Cl";
}
Cl::Cl(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 1, max_radius, Cl_initializer(), v_CM, Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0), force, Eigen::Vector3d(0.0, 0.0, 0.0), dipole = Eigen::Vector3d(0.0, 0.0, 0.0), is_linear = false, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        this->atoms[0].coor = this->coor_CM;
        this->name = "Cl";
    }
CH4::CH4() : molecule() {
    this->numatoms = 5;
    this->atoms = CH4_initializer();
    this->dipole = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->is_polyatomic = true;
    this->is_linear = false;
    this->atoms_coor = vector<Eigen::Vector3d>(this->numatoms);
    for (int i = 0; i < this->numatoms; i++) {
            this->coor_CM += this->atoms[i].coor * this->atoms[i].mass;
            this->mass += this->atoms[i].mass;
    }
    this->coor_CM /= this->mass;
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = this->atoms[i].coor - this->coor_CM;
    }
    for (int i = 0; i < numatoms; i++) {
            this->inertia_tensor(0, 0) += this->atoms[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -this->atoms[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -this->atoms[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += this->atoms[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
    }
    this->name = "CH4";
}
CH4::CH4(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, double* max_radius,\
             const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
             const Eigen::Vector3d& force, const Eigen::Vector3d& coor_CM, const Eigen::Vector3d& rotate, bool save_particles)
    : molecule(Particles, checking_lst, num_molecules, 5, max_radius, CH4_initializer(), v_CM, w_CM, torque, force, rotate, dipole = Eigen::Vector3d(0.0,0.0,0.0), is_linear = false, save_particles = save_particles) {
        this->coor_CM = coor_CM;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
        this->name = "CH4";
    }
