#ifndef MAIN_FUNC_HPP
#define MAIN_FUNC_HPP

#include <tuple>
#include <functional>
#include "specific_mol.hpp"
#include <algorithm>
#include <list>
#include <fstream>

namespace std {
    template <>
    struct hash<tuple<int, int, int>> {
        size_t operator() (const tuple<int,int,int>& t) const {
            size_t h1 = hash<int>() (get<0>(t));
            size_t h2 = hash<int>() (get<1>(t));
            size_t h3 = hash<int>() (get<2>(t));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}
namespace molcules {
    extern NO2 _NO2_;
    extern NO3 _NO3_;
    extern CO _CO_;
    extern NO _NO_;
    extern Ar _Ar_;
    extern Cl _Cl_;
    extern CH4 _CH4_;
}
void print_checking_lst(vector<vector<vector<int>>>* checking_lst);
Eigen::Vector3d E_particle(vector<molecule*>& Particles, int num_molecules, int index, Eigen::Vector3d dcoor = Eigen::Vector3d(0.0,0.0,0.0));
Eigen::Vector3d Force_particle(vector<molecule*>& Particles, int num_molecules, int index);
Eigen::Vector3d Torque_particle(vector<molecule*>& Particles, int num_molecules, int index);
template <typename moleculeType>
void map_initial_state(vector<molecule*>& Particles, int count_samples, const moleculeType& sample_molecule, int *num_molecules, double *max_radius, vector<vector<vector<int>>>* checking_lst, double Temp) {
    for (int _ = 0; _ < count_samples; _++) {
        while (true) {
            bool is_valid = true;
            double velocity = velocity_random(sample_molecule.mass, Temp);
            Eigen::Vector3d vel_ang = random_angular().normalized();
            Eigen::Vector3d angular = random_angular().normalized();
            Eigen::Vector3d omega_;
            if (!sample_molecule.is_linear) {
                omega_ = sample_molecule.get_w_CM(0.5 * sample_molecule.mass * pow(velocity_random(sample_molecule.mass, Temp), 2), angular);
            }
            else {
                omega_ = sample_molecule.get_w_CM(1.0/3.0 * sample_molecule.mass * pow(velocity_random(sample_molecule.mass, Temp), 2), angular);
            }
            new moleculeType(Particles, checking_lst, num_molecules, max_radius, velocity * vel_ang, \
                                                omega_, Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2))), random_angular(), true);
            if (Particles[Particles.size() - 1]->check_collision_wall()) {
                is_valid = false;
                Particles[Particles.size() - 1]->del(Particles, checking_lst, num_molecules);
            }
            for (int i = 0; i < *num_molecules; i++) {
                if (i != Particles[Particles.size() - 1]->index) {
                    // 개별적으로 체크
                    bool wall_collision = Particles[Particles.size() - 1]->check_collision_wall();
                    bool particle_collision = Particles[Particles.size() - 1]->check_collision(*(Particles[i]));
                    if (particle_collision || wall_collision) {
                        is_valid = false;
                        Particles[Particles.size() - 1]->del(Particles, checking_lst, num_molecules);
                        break;
                    }
                }
            }
            // print_checking_lst(checking_lst);
            if (is_valid) {
                break;
            }
        }
    }
}
template <typename Map, typename Key, typename Default>
Default getDefault(const Map& m, const Key& k, const Default& def) {
    auto it = m.find(k);
    return (it != m.end()) ? it->second : def;
}
void make_coll_list(vector<molecule*>& Particles, int num_molecules, list<vector<int>>* coll_list, double max_radius);
void coll_list_push(list<vector<int>>* coll_list, const int& index, int Particles_size);
void check_coll_wall_Particles(vector<molecule*>& Particles, const int& num_molecules);
void update_Particles(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, const int& num_molecules);
bool find_vec(const vector<int>& vec, const int& i);
/**
 * @brief check if reaction type matches, It made effective collision, etc. molecule1 + molecule2 -> molecule3 + molecule4
 * @param Particles arr of Particles(use new)
 * @param atom1_idx_vec idx of vector means effective collision atom's index
 * @param arr 0 -> reactant, 1 -> intermediate, 2 -> product (arr must ascend if not, it'll occur bug)
 * @return it reacts or not*/
template <typename molecule1, typename molecule2, typename molecule3, typename molecule4>
bool reaction12_34(vector<molecule*>& Particles, int* num_molecules, vector<vector<vector<int>>>* checking_lst, double* max_radius, list<vector<int>>* coll_list, molecule1 sample_mol1, vector<int> atom1_idx_vec, molecule2 sample_mol2, vector<int> atom2_idx_vec\
,molecule3 sample_mol3, molecule4 sample_mol4, std::array<int, 4> arr, int index1, int index2, int atom_index1, int atom_index2 , string rxn_type, string your_rxn_type, double activat_E, double delta_H, int* rxn_count, bool is_final, bool is_similar_1_3, double Temp) {
    // string particle_type[3] = {"reactant", "intermediate", "product"};
    if (rxn_type == your_rxn_type) {
        if ((Particles[index1]->get_name() == sample_mol1.get_name()) && (find_vec(atom1_idx_vec, atom_index1)) \
        && (Particles[index2]->get_name() == sample_mol2.get_name()) && (find_vec(atom2_idx_vec, atom_index2))) { //if it eff-collided
            // rxn_total++;
            if (uniform() < exp(-activat_E / (Constants::Boltzmann_Constant * Temp))) {
                // rxn_succeed++;
                if (arr[2] + arr[3] == 2) { //둘 다 중간체인 경우
                    int cnt = 0;
                    while (true) {
                        double rotate_Energy = 0.5 * sample_mol3.mass * pow(velocity_random(sample_mol3.mass, Temp), 2) * (sample_mol3.degree_of_freedom - 3) / 3;
                        Eigen::Vector3d theta = random_angular().normalized();
                        Eigen::Vector3d velocity = velocity_random(sample_mol3.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        Eigen::Vector3d angular_velocity = sample_mol3.get_w_CM(rotate_Energy, theta);
                        if (is_similar_1_3) {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        }
                        else {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index2]->coor_CM, random_angular());
                        }
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        rotate_Energy = 0.5 * sample_mol4.mass * pow(velocity_random(sample_mol4.mass, Temp), 2) * (sample_mol4.degree_of_freedom - 3) / 3;
                        theta = random_angular().normalized();
                        velocity = velocity_random(sample_mol4.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        angular_velocity = sample_mol4.get_w_CM(rotate_Energy, theta);
                        if (!is_similar_1_3) {
                            new molecule4(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        }
                        else {
                            new molecule4(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index2]->coor_CM, random_angular());
                        }
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        bool is_valid = true;
                        for (int i = 0; i < *num_molecules - 1; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        for (int i = 0; i < *num_molecules - 2; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        cnt++;
                        if (is_valid) {
                            break;
                        }
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                    }
                } else if (arr[2] + arr[3] == 3) { //molecule4는 생성물
                    int cnt = 0;
                    while (true) {
                        double rotate_Energy = 0.5 * sample_mol3.mass * pow(velocity_random(sample_mol3.mass, Temp), 2) * (sample_mol3.degree_of_freedom - 3) / 3;
                        Eigen::Vector3d theta = random_angular().normalized();
                        Eigen::Vector3d velocity = velocity_random(sample_mol3.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        Eigen::Vector3d angular_velocity = sample_mol3.get_w_CM(rotate_Energy, theta);
                        if (is_similar_1_3) {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        }
                        else {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index2]->coor_CM, random_angular());
                        }
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        bool is_valid = true;
                        for (int i = 0; i < *num_molecules - 1; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        cnt++;
                        if (is_valid) {
                            break;
                        }
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                    }
                } //둘 다 product인 경우는 굳이 고려할 필요 없이 mol1, mol2만 고려하면 됨.
                if (arr[0] + arr[1] == 0) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index1]->rotate(random_angular());
                        Particles[index1]->v_CM = Particles[index1]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index1]->w_CM = Particles[index1]->get_w_CM(Particles[index1]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index1]->check_collision_wall() || Particles[index1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    while (true) {
                        bool is_valid = true;
                        Particles[index2]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index2]->rotate(random_angular());
                        Particles[index2]->v_CM = Particles[index2]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index2]->w_CM = Particles[index2]->get_w_CM(Particles[index2]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index2) {
                                if (Particles[index2]->check_collision_wall() || Particles[index2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                } else if (arr[0] + arr[1] == 1) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index1]->rotate(random_angular());
                        Particles[index1]->v_CM = Particles[index1]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index1]->w_CM = Particles[index1]->get_w_CM(Particles[index1]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index1]->check_collision_wall() || Particles[index1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    coll_list_push(coll_list, index2, Particles.size());
                    Particles[index2]->del(Particles, checking_lst, num_molecules);
                } else {
                    if (index1 < index2) {
                        coll_list_push(coll_list, index2, Particles.size());
                        Particles[index2]->del(Particles, checking_lst, num_molecules);
                        coll_list_push(coll_list, index1, Particles.size() + 1);
                        Particles[index1]->del(Particles, checking_lst, num_molecules);
                    }
                    else {
                        coll_list_push(coll_list, index1, Particles.size());
                        Particles[index1]->del(Particles, checking_lst, num_molecules);
                        coll_list_push(coll_list, index2, Particles.size());
                        Particles[index2]->del(Particles, checking_lst, num_molecules);
                    }
                }
                if (is_final) {
                    (*rxn_count)++;
                }
                return true;
            }
            return false;
        } else if ((Particles[index2]->get_name() == sample_mol1.get_name()) && (find_vec(atom1_idx_vec, atom_index2)) \
        && (Particles[index1]->get_name() == sample_mol2.get_name()) && (find_vec(atom2_idx_vec, atom_index1))) {
            return reaction12_34(Particles, num_molecules, checking_lst, max_radius, coll_list, sample_mol1, atom1_idx_vec, sample_mol2, atom2_idx_vec, sample_mol3, sample_mol4, arr, index2, index1, atom_index2, atom_index1, rxn_type, your_rxn_type, activat_E, delta_H, rxn_count, is_final, is_similar_1_3, Temp);
        }
        return false;
    }
    return false;
}
/**
 * @brief check if reaction type matches, It made effective collision, etc. molecule1 + molecule2 -> molecule3
 * @param Particles arr of Particles(use new)
 * @param atom1_idx_vec idx of vector means effective collision atom's index
 * @param arr 0 -> reactant, 1 -> intermediate, 2 -> product (arr must ascend if not, it'll occur bug)
 * @return it reacts or not */
template <typename molecule1, typename molecule2, typename molecule3>
bool reaction12_3(vector<molecule*>& Particles, int* num_molecules, vector<vector<vector<int>>>* checking_lst, double* max_radius, list<vector<int>>* coll_list, molecule1 sample_mol1, vector<int> atom1_idx_vec, molecule2 sample_mol2, vector<int> atom2_idx_vec\
,molecule3 sample_mol3, std::array<int, 3> arr, int index1, int index2, int atom_index1, int atom_index2, string rxn_type, string your_rxn_type, double activat_E, double delta_H, int* rxn_count, bool is_final, double Temp) {
    if (rxn_type == your_rxn_type) {
        if ((Particles[index1]->get_name() == sample_mol1.get_name()) && (find_vec(atom1_idx_vec, atom_index1)) \
        && (Particles[index2]->get_name() == sample_mol2.get_name()) && (find_vec(atom2_idx_vec, atom_index2))) {
            // rxn_total++;
            if (uniform() < exp(-activat_E / (Constants::Boltzmann_Constant * Temp))) {
                // rxn_succeed++;
                if (arr[2] == 1) { //molecule3이 중간체일 경우
                    int cnt = 0;
                    while (true) {
                        double rotate_Energy = 0.5 * sample_mol3.mass * pow(velocity_random(sample_mol3.mass, Temp), 2) * (sample_mol3.degree_of_freedom - 3) / 3;
                        Eigen::Vector3d theta = random_angular().normalized();
                        Eigen::Vector3d velocity = velocity_random(sample_mol3.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        Eigen::Vector3d angular_velocity = sample_mol3.get_w_CM(rotate_Energy, theta);
                        new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), (Particles[index1]->coor_CM + Particles[index2]->coor_CM) / 2, random_angular());
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        bool is_valid = true;
                        for (int i = 0; i < *num_molecules - 1; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        cnt++;
                        if (is_valid) {
                            break;
                        }
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                    }
                } 
                if (arr[0] + arr[1] == 0) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index1]->rotate(random_angular());
                        Particles[index1]->v_CM = Particles[index1]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index1]->w_CM = Particles[index1]->get_w_CM(Particles[index1]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index1]->check_collision_wall() || Particles[index1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    while (true) {
                        bool is_valid = true;
                        Particles[index2]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index2]->rotate(random_angular());
                        Particles[index2]->v_CM = Particles[index2]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index2]->w_CM = Particles[index2]->get_w_CM(Particles[index2]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index2) {
                                if (Particles[index2]->check_collision_wall() || Particles[index2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                } else if (arr[0] + arr[1] == 1) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index1]->rotate(random_angular());
                        Particles[index1]->v_CM = Particles[index1]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index1]->w_CM = Particles[index1]->get_w_CM(Particles[index1]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index1]->check_collision_wall() || Particles[index1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    coll_list_push(coll_list, index2, Particles.size());
                    Particles[index2]->del(Particles, checking_lst, num_molecules);
                } else {
                    if (index1 < index2) {
                        coll_list_push(coll_list, index2, Particles.size());
                        Particles[index2]->del(Particles, checking_lst, num_molecules); 
                        coll_list_push(coll_list, index1, Particles.size());
                        Particles[index1]->del(Particles, checking_lst, num_molecules);
                    }
                    else {
                        coll_list_push(coll_list, index1, Particles.size());
                        Particles[index1]->del(Particles, checking_lst, num_molecules);
                        coll_list_push(coll_list, index2, Particles.size());
                        Particles[index2]->del(Particles, checking_lst, num_molecules);
                    }
                }
                if (is_final) {
                    (*rxn_count)++;
                }
                return true;
            }
            return false;
        } else if ((Particles[index2]->get_name() == sample_mol1.get_name()) && (find_vec(atom1_idx_vec, atom_index2)) \
        && (Particles[index1]->get_name() == sample_mol2.get_name()) && (find_vec(atom2_idx_vec, atom_index1))) {
            return reaction12_3(Particles, num_molecules, checking_lst, max_radius, coll_list, sample_mol1, atom1_idx_vec, sample_mol2, atom2_idx_vec, sample_mol3, arr, index2, index1, atom_index2, atom_index1, rxn_type, your_rxn_type, activat_E, delta_H, rxn_count, is_final, Temp);
        }
        return false;
    }
    return false;
}
/**
 * @brief check if reaction type matches, It made effective collision, etc. molecule1 -> molecule3 + molecule4 (But it react only when collided)
 */
template <typename molecule1, typename molecule3, typename molecule4>
bool reaction1_34(vector<molecule*>& Particles, int* num_molecules, vector<vector<vector<int>>>* checking_lst, double* max_radius, list<vector<int>>* coll_list, molecule1 sample_mol1\
,molecule3 sample_mol3, molecule4 sample_mol4, std::array<int, 3> arr, int index1, int index2, int atom_index1, int atom_index2, string rxn_type, string your_rxn_type, double activat_E, double delta_H, int* rxn_count, bool is_final, double Temp) {
    if (rxn_type == your_rxn_type) {
        if (Particles[index1]->get_name() == sample_mol1.get_name()) { //if it eff-collided
            if (uniform() < exp(-activat_E / (Constants::Boltzmann_Constant * Temp))) {
                if (arr[1] + arr[2] == 2) { //둘 다 중간체인 경우
                    int cnt = 0;
                    while (true) {
                        double rotate_Energy = 0.5 * sample_mol3.mass * pow(velocity_random(sample_mol3.mass, Temp), 2) * (sample_mol3.degree_of_freedom - 3) / 3;
                        Eigen::Vector3d theta = random_angular().normalized();
                        Eigen::Vector3d velocity = velocity_random(sample_mol3.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        Eigen::Vector3d angular_velocity = sample_mol3.get_w_CM(rotate_Energy, theta);
                        bool randomgate = uniform(Constants::get_gen()) > 0.5;
                        if (randomgate) {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        }
                        else {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index2]->coor_CM, random_angular());
                        }
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        rotate_Energy = 0.5 * sample_mol4.mass * pow(velocity_random(sample_mol4.mass, Temp), 2) * (sample_mol4.degree_of_freedom - 3) / 3;
                        theta = random_angular().normalized();
                        velocity = velocity_random(sample_mol4.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        angular_velocity = sample_mol4.get_w_CM(rotate_Energy, theta);
                        if (!randomgate) {
                            new molecule4(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        }
                        else {
                            new molecule4(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index2]->coor_CM, random_angular());
                        }
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        bool is_valid = true;
                        for (int i = 0; i < *num_molecules - 1; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        for (int i = 0; i < *num_molecules - 2; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        cnt++;
                        if (is_valid) {
                            break;
                        }
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                    }
                } else if (arr[1] + arr[2] == 3) { //molecule4는 생성물
                    int cnt = 0;
                    while (true) {
                        double rotate_Energy = 0.5 * sample_mol3.mass * pow(velocity_random(sample_mol3.mass, Temp), 2) * (sample_mol3.degree_of_freedom - 3) / 3;
                        Eigen::Vector3d theta = random_angular().normalized();
                        Eigen::Vector3d velocity = velocity_random(sample_mol3.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        Eigen::Vector3d angular_velocity = sample_mol3.get_w_CM(rotate_Energy, theta);
                        if (uniform(Constants::get_gen()) > 0.5) {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        }
                        else {
                            new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index2]->coor_CM, random_angular());
                        }
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        bool is_valid = true;
                        for (int i = 0; i < *num_molecules - 1; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        cnt++;
                        if (is_valid) {
                            break;
                        }
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                    }
                } //둘 다 product인 경우는 굳이 고려할 필요 없이 mol1만 고려하면 됨.
                if (arr[0] == 0) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index1]->rotate(random_angular());
                        Particles[index1]->v_CM = Particles[index1]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index1]->w_CM = Particles[index1]->get_w_CM(Particles[index1]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index1]->check_collision_wall() || Particles[index1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    while (true) {
                        bool is_valid = true;
                        Particles[index2]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index2]->rotate(random_angular());
                        Particles[index2]->v_CM = Particles[index2]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index2]->w_CM = Particles[index2]->get_w_CM(Particles[index2]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index2) {
                                if (Particles[index2]->check_collision_wall() || Particles[index2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                } else if (arr[0] == 1) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index2]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index2]->rotate(random_angular());
                        Particles[index2]->v_CM = Particles[index2]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index2]->w_CM = Particles[index2]->get_w_CM(Particles[index2]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index2]->check_collision_wall() || Particles[index2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    coll_list_push(coll_list, index1, Particles.size());
                    Particles[index1]->del(Particles, checking_lst, num_molecules);
                }
                if (is_final) {
                    (*rxn_count)++;
                }
                return true;
            }
            return false;
        } else if ((Particles[index2]->get_name() == sample_mol1.get_name())) {
            return reaction1_34(Particles, num_molecules, checking_lst, max_radius, coll_list, sample_mol1, sample_mol3, sample_mol4, arr, index2, index1, atom_index2, atom_index1, rxn_type, your_rxn_type, activat_E, delta_H, rxn_count, is_final, Temp);
        }
        return false;
    }
    return false;
}
/**
 *@brief molecule1 -> molecule3 (But it reat only when collided)  
 */
template <typename molecule1, typename molecule3>
bool reaction1_3(vector<molecule*>& Particles, int* num_molecules, vector<vector<vector<int>>>* checking_lst, double* max_radius, list<vector<int>>* coll_list, molecule1 sample_mol1\
,molecule3 sample_mol3, std::array<int, 2> arr, int index1, int index2, int atom_index1, int atom_index2, string rxn_type, string your_rxn_type, double activat_E, double delta_H, int* rxn_count, bool is_final, double Temp) {
    if (rxn_type == your_rxn_type) {
        if (Particles[index1]->get_name() == sample_mol1.get_name()) {
            if (uniform() < exp(-activat_E / (Constants::Boltzmann_Constant * Temp))) {
                if (arr[1] == 1) {
                    int cnt = 0;
                    while (true) {
                        double rotate_Energy = 0.5 * sample_mol3.mass * pow(velocity_random(sample_mol3.mass, Temp), 2) * (sample_mol3.degree_of_freedom - 3) / 3;
                        Eigen::Vector3d theta = random_angular().normalized();
                        Eigen::Vector3d velocity = velocity_random(sample_mol3.mass, Temp) * theta;
                        theta = random_angular().normalized();
                        Eigen::Vector3d angular_velocity = sample_mol3.get_w_CM(rotate_Energy, theta);
                        new molecule3(Particles, checking_lst, num_molecules, max_radius, velocity, angular_velocity, Eigen::Vector3d(0.0,0.0,0.0), Eigen::Vector3d(0.0,0.0,0.0), Particles[index1]->coor_CM, random_angular());
                        if (cnt >= 100) {
                            Particles[*num_molecules - 1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        }
                        bool is_valid = true;
                        for (int i = 0; i < *num_molecules - 1; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[*num_molecules - 1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        cnt++;
                        if (is_valid) {
                            break;
                        }
                        Particles[*num_molecules - 1]->del(Particles, checking_lst, num_molecules);
                    }
                }
                if (arr[0] == 0) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index1]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index1]->rotate(random_angular());
                        Particles[index1]->v_CM = Particles[index1]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index1]->w_CM = Particles[index1]->get_w_CM(Particles[index1]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index1]->check_collision_wall() || Particles[index1]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    while (true) {
                        bool is_valid = true;
                        Particles[index2]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index2]->rotate(random_angular());
                        Particles[index2]->v_CM = Particles[index2]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index2]->w_CM = Particles[index2]->get_w_CM(Particles[index2]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index2) {
                                if (Particles[index2]->check_collision_wall() || Particles[index2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                } else if (arr[0] == 1) {
                    while (true) {
                        bool is_valid = true;
                        Particles[index2]->coor_CM = Eigen::Vector3d(uniform(Constants::get_gen(), 0.0, Constants::SCREEN(0)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(1)), uniform(Constants::get_gen(), 0.0, Constants::SCREEN(2)));
                        Particles[index2]->rotate(random_angular());
                        Particles[index2]->v_CM = Particles[index2]->v_CM.norm() * random_angular().normalized();
                        Eigen::Vector3d theta = random_angular().normalized();
                        Particles[index2]->w_CM = Particles[index2]->get_w_CM(Particles[index2]->Rotation_Energy(), theta);
                        for (int i = 0; i < *num_molecules; i++) {
                            if (i != index1 && i != index2) {
                                if (Particles[index2]->check_collision_wall() || Particles[index2]->check_collision(*(Particles[i]))) {
                                    is_valid = false;
                                    break;
                                }
                            }
                        }
                        if (is_valid) {
                            break;
                        }
                    }
                    coll_list_push(coll_list, index1, Particles.size());
                    Particles[index1]->del(Particles, checking_lst, num_molecules);
                }
                if (is_final) {
                    (*rxn_count)++;
                }
                return true;
            }
            return false;
        } else if ((Particles[index2]->get_name() == sample_mol1.get_name())) {
            return reaction1_3(Particles, num_molecules, checking_lst, max_radius, coll_list, sample_mol1, sample_mol3, arr, index2, index1, atom_index2, atom_index1, rxn_type, your_rxn_type, activat_E, delta_H, rxn_count, is_final, Temp);
        }
        return false;
    }
    return false;
}
void initialize(vector<molecule*>& Particles, int* num_molecules, vector<vector<vector<int>>>* checking_lst, double* max_radius, list<vector<int>>* coll_list);
vector<int> pop_vec(list<vector<int>>* coll_list);
void average_Kinetic(vector<molecule*> Particles);
void print_coll_list(list<vector<int>>* coll_list);
void save_energy_dist(vector<molecule*> Particles);
#endif // MAIN_FUNC_HPP