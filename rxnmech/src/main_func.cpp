#include "main_func.hpp"
namespace molcules {
    NO2 _NO2_ = NO2();
    NO3 _NO3_ = NO3();
    CO _CO_ = CO();
    NO _NO_ = NO();
    Ar _Ar_ = Ar();
    Cl _Cl_ = Cl();
    CH4 _CH4_ = CH4();
}
void print_checking_lst(vector<vector<vector<int>>>* checking_lst) {
    cout << "Checking list: " << endl;
    for (size_t i = 0; i < (*checking_lst).size(); i++) {
        for (size_t j = 0; j < (*checking_lst)[i].size(); j++) {
            cout << "(" << (*checking_lst)[i][j][0] << ", " << (*checking_lst)[i][j][1] << ")";
        }
        cout << endl;
    }
}
Eigen::Vector3d E_particle(vector<molecule*>& Particles, int num_molecules, int index, Eigen::Vector3d dcoor) {
    Eigen::Vector3d E = Eigen::Vector3d(0.0, 0.0, 0.0);
    for (size_t i = 0; i < Particles.size(); i++) {
        if (static_cast<int>(i) == index) continue;
        E += E_field(Particles[i]->dipole, Particles[index]->coor_CM + dcoor - Particles[i]->coor_CM, Constants::Coulomb_Constant);
    }
    return E;
}
Eigen::Vector3d Force_particle(vector<molecule*>& Particles, int num_molecules, int index) {
    double dx = Constants::dcoor;
    double dy = Constants::dcoor;
    double dz = Constants::dcoor;
    double dUdx = (Potential_E(Particles[index]->dipole, E_particle(Particles, num_molecules, index, Eigen::Vector3d(dx, 0.0, 0.0))) - 
                   Potential_E(Particles[index]->dipole, E_particle(Particles, num_molecules, index))) / dx;
    double dUdy = (Potential_E(Particles[index]->dipole, E_particle(Particles, num_molecules, index, Eigen::Vector3d(0.0, dy, 0.0))) -
                   Potential_E(Particles[index]->dipole, E_particle(Particles, num_molecules, index))) / dy;
    double dUdz = (Potential_E(Particles[index]->dipole, E_particle(Particles, num_molecules, index, Eigen::Vector3d(0.0, 0.0, dz))) -
                   Potential_E(Particles[index]->dipole, E_particle(Particles, num_molecules, index))) / dz;
    return Eigen::Vector3d(-dUdx, -dUdy, -dUdz);
}
Eigen::Vector3d Torque_particle(vector<molecule*>& Particles, int num_molecules, int index) {
    return Torque(Particles[index]->dipole, E_particle(Particles, num_molecules, index));
}
void make_coll_list(vector<molecule*>& Particles, int num_molecules, list<vector<int>>* coll_list, double max_radius) {
    (*coll_list).clear();
    unordered_map<tuple<int, int, int>, vector<int>> gridmap = {};
    tuple<int, int, int> grid;
    vector<int> gridElements;
    vector<int> gridSurrElements;
    for (int i = 0; i < static_cast<int>(Particles.size()); i++) {
        grid = tuple<int, int, int>((int)(Particles[i]->coor_CM(0)/(2 * max_radius)), (int)(Particles[i]->coor_CM(1)/(2 * max_radius)), (int)(Particles[i]->coor_CM(2)/(2 * max_radius)));
        gridElements = getDefault(gridmap, grid, vector<int>{});
        gridElements.push_back(i);
        gridmap[grid] = gridElements;
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                for (int l = -1; l < 2; l++) {
                    gridSurrElements = getDefault(gridmap, tuple<int, int, int>(get<0>(grid) + j, get<1>(grid) + k, get<2>(grid) + l), vector<int>{});
                    for (int m = 0; m < static_cast<int>(gridSurrElements.size()); m++) {
                        if (i != gridSurrElements[m]) {
                            (*coll_list).push_back(vector<int>{i, gridSurrElements[m]});
                        }
                    }
                }
            }
        }
    }
}
// map_initial_state(Particles, 20, molcules::_NO2_, num_molecules, max_radius, checking_lst);
void coll_list_push(list<vector<int>>* coll_list, const int& index, int Particles_size) {
    auto it = coll_list->begin();
    while (it != coll_list->end()) {
        bool erase_flag = false;
        if ((*it)[0] == index || (*it)[1] == index) {
            erase_flag = true;
        } else if ((*it)[0] == Particles_size - 1) {
            (*it)[0] = index;
        } else if ((*it)[1] == Particles_size - 1) {
            (*it)[1] = index;
        }
        if (erase_flag) {
            it = coll_list->erase(it);  // erase 후 iterator 반환됨
        } else {
            ++it;
        }
    }
}
void check_coll_wall_Particles(vector<molecule*>& Particles, const int& num_molecules) {
    for (size_t i = 0; i < Particles.size(); i++) {
        Particles[i]->collision_wall();
    }
}
void update_Particles(vector<molecule*>& Particles, vector<vector<vector<int>>> *checking_lst, const int& num_molecules) {
    for (size_t i = 0; i < Particles.size(); i++) {
        Particles[i]->update(checking_lst, num_molecules, Constants::dt);
    }
}
bool find_vec(const vector<int>& vec, const int& i) {
    for (size_t j = 0; j < vec.size(); j++) {
        if (vec[j] == i) {
            return true;
        }
    }
    return false;
}
void initialize(vector<molecule*>& Particles, int* num_molecules, vector<vector<vector<int>>>* checking_lst, double* max_radius, list<vector<int>>* coll_list) {
    int len = (int) Particles.size();
    for (int i = 0; i < len; i++) {
        Particles[Particles.size() - 1]->del(Particles, checking_lst, num_molecules);
    }
    (*checking_lst).clear();
    *max_radius = 0.0;
    (*coll_list).clear();
}
vector<int> pop_vec(list<vector<int>>* coll_list) {
    if ((*coll_list).empty()) {
        return vector<int>();
    } else {
        vector<int> popped = (*coll_list).front();
        (*coll_list).pop_front();
        return popped;
    }
}
void average_Kinetic(vector<molecule*> Particles) {
    double total_KE = 0.0;
    for (size_t i = 0; i < Particles.size(); i++) {
        total_KE += Particles[i]->Kinetic_Energy();
    }
    double average_KE = total_KE / Particles.size();
    cout << "Average Kinetic Energy: " << average_KE << endl;
}
void print_coll_list(list<vector<int>>* coll_list) {
    cout << "Collision List: " << endl;
    for (const auto& pair : *coll_list) {
        cout << "(" << pair[0] << ", " << pair[1] << ")" << ", ";
    }
    cout << endl;
}
void save_energy_dist(vector<molecule*> Particles) {
    ofstream file("../energy_distribution.txt", std::ios::app);
    if (!file.is_open()) {
        cerr << "Error opening file for writing energy distribution." << endl;
        return;
    }
    for (const auto& particle : Particles) {
        file << particle->Kinetic_Energy() << " ";
    }
    file << endl;
    file.close();
}