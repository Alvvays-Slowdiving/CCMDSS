#include "main_func.hpp"
#include <numeric>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>
int main_NO2(int NO2__, int CO__, int Ar__, double Temp__, int rxn____, int number__);
int main_CH4(int CH4__, int Cl__, double Temp__, int rxn____, int number__);

int main() {
    clock_t start, end;
    double duration;
    start = clock();
    for (double T = 10000.0; T <= 20000.0; T += 1000.0) {
        for (int i = 30; i <= 90; i += 5) {
            main_NO2(i, 30, 0, T, 1000, 10);
            cout << "----------------------------------------" << endl;
        }
    }
    for (int i = 30; i <= 90; i += 5) {
        main_NO2(50, i, 20, 15000.0, 1000, 10);
        cout << "----------------------------------------" << endl;
    }
    // for (double T = 10000.0; T <= 20000.0; T += 1000.0) {
    //     for (int i = 30; i <= 90; i += 5) {
    //         main_CH4(i, 50, T, 1000, 10);
    //     }
    //     for (int i = 30; i <= 90; i += 5) {
    //         main_CH4(50, i, T, 1000, 10);
    //     }
    // }  
    end = clock();
    duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Total execution time: " << duration << " seconds" << endl;
    return 0;
}
int main_CH4(int CH4__, int Cl__, double Temp__, int rxn____, int number__) {
    int CH4_ = CH4__;
    int Cl_ = Cl__;
    double Temp = Temp__;
    int rxn__ = rxn____;
    int number = number__;
    string rxn_type = "CH4 + Cl -> CH3 + HCl";

    vector<molecule*> Particles = vector<molecule*>();
    int number_of_particles = 0;
    int* num_molecules = &number_of_particles;
    vector<vector<vector<int>>> init_checking_lst = vector<vector<vector<int>>> {};
    vector<vector<vector<int>>>* checking_lst = &init_checking_lst;
    double initial_radius = 0.0;
    double* max_radius = &initial_radius;
    list<vector<int>> init_coll_list = list<vector<int>> {{}};
    list<vector<int>>* coll_list = &init_coll_list;
    vector<double> rxn_rate_lst;
    double t;
    for (int i = 0; i < number; i++) {
        int cnt = 0;
        map_initial_state(Particles, CH4_, molcules::_CH4_, num_molecules, max_radius, checking_lst, Temp);
        map_initial_state(Particles, Cl_, molcules::_Cl_, num_molecules, max_radius, checking_lst, Temp);
        t = 0;
        int init_rxn_count = 0;
        int* rxn_count = &init_rxn_count;
        double rxn_rate = 0;

        while (true) {
            check_coll_wall_Particles(Particles, *num_molecules);
            make_coll_list(Particles, *num_molecules, coll_list, *max_radius);
            // list<vector<int>> copy_coll_list = *coll_list;
            while (true) {
                vector<int> coll_pair = pop_vec(coll_list);
                if (coll_pair.empty()) {
                    break;
                }
                int index1 = coll_pair[0];
                int index2 = coll_pair[1];
                int atom_index1 = -1;
                int atom_index2 = -1;
                bool is_valid = Particles[index1]->collision(checking_lst, *(Particles[index2]), &atom_index1, &atom_index2);
                if (is_valid) {
                    if (reaction12_34(Particles, num_molecules, checking_lst, max_radius, coll_list, molcules::_CH4_, vector<int>{1, 2, 3, 4}, molcules::_Cl_, vector<int>{0}, molcules::_Ar_, molcules::_Ar_, {0, 0, 2, 2}, index1, index2, atom_index1, atom_index2, rxn_type, "CH4 + Cl -> CH3 + HCl", Constants::Ea_5, 0.0, rxn_count, true, false, Temp)) {
                        double percentage = (static_cast<double>(*rxn_count) / rxn__) * 100.0;
                        int bar_width = 50; // 프로그레스 바의 전체 너비
                        int filled_width = static_cast<int>(std::floor(bar_width * percentage / 100.0));
                        std::string bar_filled(filled_width, '#'); // 채워진 부분
                        std::string bar_empty(bar_width - filled_width, '-'); // 비어있는 부분
                        cout << fixed << setprecision(2) << percentage << "% "
                             << "[" << bar_filled << bar_empty << "] "
                             << *rxn_count << "/" << rxn__ << "\r" << flush;
                        
                        continue;
                    }
                }
            }
            update_Particles(Particles, checking_lst, *num_molecules);
            t += Constants::dt;

            if (*rxn_count >= rxn__) {
                cout << string(80, ' ') << "\r";
                rxn_rate = static_cast<double>(*rxn_count) / t;
                rxn_rate_lst.push_back(rxn_rate / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0));
                cout << fixed << setprecision(10);
                cout << "rxn rate: " << rxn_rate / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0) << "M/s" << endl;
                initialize(Particles, num_molecules,checking_lst, max_radius, coll_list);
                break;
            }
        }
    }
    cout << fixed << setprecision(10);
    cout << "CH4 concentration: " << static_cast<double>(CH4_) / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0) << "M" << endl;
    cout << "Cl concentration: " << static_cast<double>(Cl_) / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0) << "M" << endl;
    cout << "Temperature: " << Temp << "K" << endl;
    double mean_rate = accumulate(rxn_rate_lst.begin(), rxn_rate_lst.end(), 0.0) / static_cast<double>(rxn_rate_lst.size());
    cout << "mean reaction rate: " << mean_rate << " M/s" << endl;
    double std_dev = sqrt(inner_product(rxn_rate_lst.begin(), rxn_rate_lst.end(), rxn_rate_lst.begin(), 0.0) / static_cast<double>(rxn_rate_lst.size()) - pow(mean_rate, 2));
    cout << "standard error: " << std_dev / sqrt(static_cast<double>(rxn_rate_lst.size())) << " M/s"<<endl;
    return 0; 
}
int main_NO2(int NO2__, int CO__, int Ar__, double Temp__, int rxn____, int number__) {
    // 기본값 설정
    int NO2_ = NO2__;
    int CO_ = CO__;
    int Ar_ = Ar__;
    double Temp = Temp__;
    int rxn__ = rxn____;
    int number = number__;
    string rxn_type = "NO2 + CO -> NO + CO2"; 
    // if (argc > 1) NO2_ = atoi(argv[1]);
    // if (argc > 2) CO_ = atoi(argv[2]);
    // if (argc > 3) Ar_ = atoi(argv[3]);
    // if (argc > 4) Temp = atof(argv[4]);
    // if (argc > 5) rxn__ = atoi(argv[5]);
    // if (argc > 6) number = atoi(argv[6]);
    
    vector<molecule*> Particles = vector<molecule*>();
    int number_of_particles = 0;
    int* num_molecules = &number_of_particles;
    vector<vector<vector<int>>> init_checking_lst = vector<vector<vector<int>>> {};
    vector<vector<vector<int>>>* checking_lst = &init_checking_lst;
    double initial_radius = 0.0;
    double* max_radius = &initial_radius;
    list<vector<int>> init_coll_list = list<vector<int>> {{}};
    list<vector<int>>* coll_list = &init_coll_list;
    vector<double> rxn_rate_lst;
    double t;
    for (int i = 0; i < number; i++) {
        int cnt = 0;
        map_initial_state(Particles, NO2_, molcules::_NO2_, num_molecules, max_radius, checking_lst, Temp);
        map_initial_state(Particles, CO_, molcules::_CO_, num_molecules, max_radius, checking_lst, Temp);
        map_initial_state(Particles, Ar_, molcules::_Ar_, num_molecules, max_radius, checking_lst, Temp);
        t = 0;
        int init_rxn_count = 0;
        int* rxn_count = &init_rxn_count;
        double rxn_rate = 0;
        while (true) {
            check_coll_wall_Particles(Particles, *num_molecules);
            make_coll_list(Particles, *num_molecules, coll_list, *max_radius);
            // list<vector<int>> copy_coll_list = *coll_list;
            while (true) {
                vector<int> coll_pair = pop_vec(coll_list);
                if (coll_pair.empty()) {
                    break;
                }
                int index1 = coll_pair[0];
                int index2 = coll_pair[1];
                int atom_index1 = -1;
                int atom_index2 = -1;
                bool is_valid = Particles[index1]->collision(checking_lst, *(Particles[index2]), &atom_index1, &atom_index2);
                if (is_valid) {
                    if (reaction12_34(Particles, num_molecules, checking_lst, max_radius, coll_list, molcules::_NO2_, vector<int>{0}, molcules::_NO2_, vector<int>{1,2}, molcules::_NO3_, molcules::_NO_, {0, 0, 1, 2}, index1, index2, atom_index1, atom_index2, rxn_type, "NO2 + CO -> NO + CO2", Constants::Ea_1, 0.0, rxn_count, false, true, Temp)) {
                        continue;
                    }
                    if (reaction12_3(Particles, num_molecules, checking_lst, max_radius, coll_list, molcules::_CO_, vector<int>{0}, molcules::_NO3_, vector<int>{1, 2, 3},  molcules::_NO2_, {0, 1, 2}, index1, index2, atom_index1, atom_index2, rxn_type, "NO2 + CO -> NO + CO2", Constants::Ea_3, 0.0, rxn_count, true, Temp)) {
                        double percentage = (static_cast<double>(*rxn_count) / rxn__) * 100.0;
                        int bar_width = 50; // 프로그레스 바의 전체 너비
                        int filled_width = static_cast<int>(std::floor(bar_width * percentage / 100.0));
                        std::string bar_filled(filled_width, '#'); // 채워진 부분
                        std::string bar_empty(bar_width - filled_width, '-'); // 비어있는 부분
                        cout << fixed << setprecision(2) << percentage << "% "
                             << "[" << bar_filled << bar_empty << "] "
                             << *rxn_count << "/" << rxn__ << "\r" << flush;
                        
                        continue;
                    }
                }
            }
            update_Particles(Particles, checking_lst, *num_molecules);
            t += Constants::dt;
            // if (t > 0.1e-10) {
            //     cnt++;
            //     if (cnt >= 100) {
            //         exit(0);
            //     }
            //     // cout << cnt <<", Particles: " << Particles.size() << ", checking_lst: " << checking_lst->size() << endl;
            //     // average_Kinetic(Particles);
            //     t -= 0.1e-10;
            // }
            if (*rxn_count >= rxn__) {
                cout << string(80, ' ') << "\r" << flush;
                rxn_rate = static_cast<double>(*rxn_count) / t;
                rxn_rate_lst.push_back(rxn_rate / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0));
                cout << fixed << setprecision(10);
                cout << "rxn rate: " << rxn_rate / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0) << "M/s" << endl;
                initialize(Particles, num_molecules,checking_lst, max_radius, coll_list);
                break;
            }
        }
    }
    cout << fixed << setprecision(10);
    cout << "NO2 concentration: " << static_cast<double>(NO2_) / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0) << "M" << endl;
    cout << "CO concentration: " << static_cast<double>(CO_) / (Constants::Avogadro_Number * Constants::SCREEN(0) * Constants::SCREEN(1) * Constants::SCREEN(2) * 1000.0) << "M" << endl;
    cout << "Temperature: " << Temp << "K" << endl;
    double mean_rate = accumulate(rxn_rate_lst.begin(), rxn_rate_lst.end(), 0.0) / static_cast<double>(rxn_rate_lst.size());
    cout << "mean reaction rate: " << mean_rate << " M/s" << endl;
    double std_dev = sqrt(inner_product(rxn_rate_lst.begin(), rxn_rate_lst.end(), rxn_rate_lst.begin(), 0.0) / static_cast<double>(rxn_rate_lst.size()) - pow(mean_rate, 2));
    cout << "standard error: " << std_dev / sqrt(static_cast<double>(rxn_rate_lst.size())) << " M/s"<<endl;
    return 0;
}