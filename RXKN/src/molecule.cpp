#include "molecule.hpp"
molecule::molecule()
{
    this->dipole = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->atoms = vector<atom>{};
    this->numatoms = 0;
    this->v_CM = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->coor_CM = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->w_CM = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->torque = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->force = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->atoms_coor = vector<Eigen::Vector3d>{};
    this->inertia_tensor = Eigen::Matrix3d::Zero();
    this->mass = 0.0;
    this->is_polyatomic = false;
    this->is_linear = false;
    this->linear_axis = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->degree_of_freedom = 6; // Default for a molecule
    this->radius = -1.0; // Default radius
    this->index = -1; // Default index
    this->name = "";
}
molecule::molecule(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules, int num, double* max_radius, vector<atom> atoms_array,\
                   const Eigen::Vector3d& v_CM, const Eigen::Vector3d& w_CM, const Eigen::Vector3d& torque,\
                   const Eigen::Vector3d& force, const Eigen::Vector3d& rotate, const Eigen::Vector3d& dipole,\
                   bool is_linear, bool save_particles)
    {
        this->dipole = dipole;
        this->atoms.reserve(num);
        this->atoms.clear();
        this->atoms_coor.reserve(num);
        this->atoms_coor.clear();
        for (int i = 0; i < num; i++) {
            this->atoms.push_back(atoms_array[i]);
        }
        this->coor_CM = Eigen::Vector3d(0.0, 0.0, 0.0);
        this->mass = 0.0;
        this->numatoms = num;
        this->v_CM = v_CM;
        for (int i = 0; i < num; i++) {
            this->coor_CM += atoms_array[i].coor * atoms_array[i].mass;
            this->mass += atoms_array[i].mass;
        }
        this->coor_CM /= this->mass;
        this->w_CM = w_CM;
        this->torque = torque;
        this->force = force;
        for (int i = 0; i < num; i++) {
            this->atoms_coor.push_back(atoms_array[i].coor - this->coor_CM);
        }
        this->inertia_tensor = Eigen::Matrix3d::Zero();
        for (int i = 0; i < num; i++) {
            this->inertia_tensor(0, 0) += atoms_array[i].mass * (this->atoms_coor[i](1) * this->atoms_coor[i](1) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(0, 1) += -atoms_array[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(0, 2) += -atoms_array[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(1, 0) += -atoms_array[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](1);
            this->inertia_tensor(1, 1) += atoms_array[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](2) * this->atoms_coor[i](2));
            this->inertia_tensor(1, 2) += -atoms_array[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 0) += -atoms_array[i].mass * this->atoms_coor[i](0) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 1) += -atoms_array[i].mass * this->atoms_coor[i](1) * this->atoms_coor[i](2);
            this->inertia_tensor(2, 2) += atoms_array[i].mass * (this->atoms_coor[i](0) * this->atoms_coor[i](0) + this->atoms_coor[i](1) * this->atoms_coor[i](1));
        }
        this->is_polyatomic = (bool) (num > 1);
        this->is_linear = is_linear;
        if (is_linear) {
            for (int i = 0; i < num; i++) {
                if (!this->atoms_coor[i].isZero(1e-12)) {  // 1e-12 허용오차로 비교
                    this->linear_axis = this->atoms_coor[i].normalized();  // 방향 벡터로 정규화
                    break;
                }
            }
        }
        if (this->is_polyatomic) {
            if (this->is_linear) {
                this->degree_of_freedom = 5; // Linear molecules have 5 degrees of freedom
            } else {
                this->degree_of_freedom = 6; // Non-linear molecules have 6 degrees of freedom
            }
        }
        else {
            this->degree_of_freedom = 3; // Monoatomic molecules have 3 degrees of freedom
        }
        Eigen::Matrix3d Rotat = Rotator_w(rotate, 1);
        if (this->is_linear) {
            this->linear_axis = Rotat * this->linear_axis;  // 회전 후 선형 축 업데이트
        }
        this->inertia_tensor = Rotat * this->inertia_tensor * Rotat.transpose();
        this->dipole = Rotat * this->dipole;  // 회전 후 쌍극자 모멘트 업데이트
        this->w_CM = Rotat * this->w_CM;  // 회전 후 각속도 업데이트
        this->v_CM = Rotat * this->v_CM;  // 회전 후 중심 속도 업데이트
        double max_radius_value = 0.0;
        for (int i = 0; i < num; i++) {
            this->atoms_coor[i] = Rotat * this->atoms_coor[i];
            if (this->atoms_coor[i].norm() + this->atoms[i].radius > max_radius_value) {
                max_radius_value = this->atoms_coor[i].norm() + this->atoms[i].radius;
            }
        }
        this->radius = max_radius_value;
        if (*max_radius < max_radius_value) {
            *max_radius = max_radius_value;
        }
        
        if (save_particles) {
            this->index = *num_molecules;
            Particles.push_back(this);
            *num_molecules = Particles.size();
            (*checking_lst).push_back(vector<vector<int>>{});
            for (size_t i = 0; i < Particles.size() - 1; i++) {
                (*checking_lst)[i].push_back(vector<int>{-1, -1});
            }
        } else {
            this->index = -1;
        }
        this->name = "";
    }
double molecule::Kinetic_Energy() const {
    return 0.5 * this->mass * this->v_CM.squaredNorm() + 0.5 * this->w_CM.dot(this->inertia_tensor * this->w_CM);
}
double molecule::trans_Energy() const {
    return 0.5 * this->mass * this->v_CM.squaredNorm();
}
double molecule::Rotation_Energy() const {
    return 0.5 * this->w_CM.dot(this->inertia_tensor * this->w_CM);
}
void molecule::rotate(Eigen::Vector3d omega) {
    Eigen::Matrix3d Rotat = Rotator_w(omega, 1);
    if (this->is_linear) {
        this->linear_axis = Rotat * this->linear_axis;  // 회전 후 선형 축 업데이트
    }
    this->inertia_tensor = Rotat * this->inertia_tensor * Rotat.transpose();
    this->dipole = Rotat * this->dipole;  // 회전 후 쌍극자 모멘트 업데이트
    this->w_CM = Rotat * this->w_CM;  // 회전 후
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = Rotat * this->atoms_coor[i];
    }
}
void molecule::rotation(double dt = Constants::dt) {
    Eigen::Matrix3d Rotat = Rotator_w(this->w_CM, dt);
    if (this->is_linear) {
        this->linear_axis = Rotat * this->linear_axis;  // 회전 후 선형 축 업데이트
    }
    this->inertia_tensor = Rotat * this->inertia_tensor * Rotat.transpose();
    this->dipole = Rotat * this->dipole;  // 회전 후 쌍극자 모멘트 업데이트
    for (int i = 0; i < this->numatoms; i++) {
        this->atoms_coor[i] = Rotat * this->atoms_coor[i];
        this->atoms[i].coor = Rotat * this->atoms[i].coor;
    }
}
void molecule::update(vector<vector<vector<int>>>* checking_lst, int num_particles, double dt = Constants::dt) {
    bool boolean = false; // is_collided
    for (int i = 0; i < this->index; i++ ) {
        if ((*checking_lst)[i][this->index - i - 1][0] != -1) {
            boolean = true;
        }
    }
    for (int i = this->index + 1; i < num_particles; i++) {
        if ((*checking_lst)[this->index][i - this->index - 1][0] != -1) {
            boolean = true;
        }
    }
    if (!boolean) {
        if (this->is_polyatomic && !this->is_linear) {
            Eigen::Vector3d ang_acc = this->inertia_tensor.inverse() * (this->torque); /* -this->w_CM.cross(this->inertia_tensor * this->w_CM) */
            this->w_CM += ang_acc * dt;
            Eigen::Vector3d acc = this->force / this->mass;
            this->v_CM += acc * dt;
        }
        else if (this->is_polyatomic) {
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(this->linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(this->linear_axis) * this->linear_axis).normalized();
            double I_CM = crossVector.dot(this->inertia_tensor * crossVector);
            this->w_CM += (this->torque / I_CM) * dt;
            Eigen::Vector3d acc = this->force / this->mass;
            this->v_CM += acc * dt;
        }
        else {
            Eigen::Vector3d acc = this->force / this->mass;
            this->v_CM += acc * dt;
        }
        this->rotation();
        this->coor_CM += this->v_CM * dt;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
    }
    else {
        this->rotation();
        this->coor_CM += this->v_CM * dt;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
    }
}
bool molecule::check_collision(const molecule& other) const {
    bool boolean = false;
    if (this->radius + other.radius >= (this->coor_CM - other.coor_CM).norm()) {
        for (int i = 0; i < this->numatoms; i++) {
            for (int j = 0; j < other.numatoms; j++) {
                if (this->atoms[i].checkCollision(other.atoms[j])) {
                    boolean = true;
                    break;
                }
            }
            if (boolean) break;
        }
    }
    return boolean;
}
bool molecule::collision(vector<vector<vector<int>>>* checking_lst, molecule& other, int* i, int* j) {
    bool boolean = false;
    if ((this->radius + other.radius) >= (this->coor_CM - other.coor_CM).norm()) {
        for (int i1 = 0; i1 < this->numatoms; i1++) {
            for (int j1 = 0; j1 < other.numatoms; j1++) {
                if (this->index < other.index) {
                    if ((*checking_lst)[this->index][other.index - this->index - 1] != vector<int>{i1, j1}) {
                        if (this->atoms[i1].checkCollision(other.atoms[j1])) {
                            (*checking_lst)[this->index][other.index - this->index - 1] = vector<int>{i1, j1};
                            boolean = true;
                            this->cal_collision(other, i1, j1);
                            *i = i1;
                            *j = j1;
                            break;
                        }
                    }
                }
                else {
                    if ((*checking_lst)[other.index][this->index - other.index - 1] != vector<int>{j1, i1}) {
                        if (this->atoms[i1].checkCollision(other.atoms[j1])) {
                            (*checking_lst)[other.index][this->index - other.index - 1] = vector<int>{j1, i1};
                            boolean = true;
                            this->cal_collision(other, i1, j1);
                            *i = i1;
                            *j = j1;
                            break;
                        }
                    }
                }
            }
            if (boolean) {
                break;
            }
        }
    }
    if (!boolean) {
        *i = -1;
        *j = -1;
        if (this->index < other.index) {
            (*checking_lst)[this->index][other.index - this->index - 1] = vector<int>{-1, -1};
        }
        else {
            (*checking_lst)[other.index][this->index - other.index - 1] = vector<int>{-1, -1};
        }
    }
    return boolean;
}
void molecule::cal_collision(molecule& other, const int& i, const int& j) {
    if (this->is_polyatomic && other.is_polyatomic) {
        if ((!this->is_linear) && (!other.is_linear)) {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r1 = this->atoms[i].coor - this->coor_CM;
            Eigen::Vector3d r2 = other.atoms[j].coor - other.coor_CM;
            double m_inv = 1/this->mass + 1/other.mass + (r1.cross(N_hat)).transpose() * this->inertia_tensor.inverse() * (r1.cross(N_hat)) + r2.cross(N_hat).transpose() * other.inertia_tensor.inverse() * (r2.cross(N_hat));
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) + 2 * r2.cross(N_hat).dot(other.w_CM) - 2 * r1.cross(N_hat).dot(this->w_CM);
            double Impulse = vel / m_inv;
            this->v_CM +=  Impulse * N_hat / this->mass;
            this->w_CM += Impulse * (this->inertia_tensor.inverse() * r1.cross(N_hat));
            other.v_CM -= Impulse * N_hat / other.mass;
            other.w_CM -= Impulse * (other.inertia_tensor.inverse() * r2.cross(N_hat));
        }
        else if (!other.is_linear) {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r1 = this->atoms[i].coor - this->coor_CM;
            Eigen::Vector3d r2 = other.atoms[j].coor - other.coor_CM;
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(this->linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(this->linear_axis) * this->linear_axis).normalized();
            double I_CM = crossVector.dot(this->inertia_tensor * crossVector);
            double m_inv = 1/this->mass + 1/other.mass + (r1.cross(N_hat).squaredNorm() / I_CM) + r2.cross(N_hat).transpose() * other.inertia_tensor.inverse() * (r2.cross(N_hat));
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) + 2 * r2.cross(N_hat).dot(other.w_CM) - 2 * this->w_CM.dot(r1.cross(N_hat));
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * r1.cross(N_hat) / I_CM;
            other.v_CM -= Impulse * N_hat / other.mass;
            other.w_CM -= Impulse * (other.inertia_tensor.inverse() * r2.cross(N_hat));
        }
        else if (!this->is_linear) {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r1 = this->atoms[i].coor - this->coor_CM;
            Eigen::Vector3d r2 = other.atoms[j].coor - other.coor_CM;
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(other.linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(other.linear_axis) * other.linear_axis).normalized();
            double I_CM = crossVector.dot(other.inertia_tensor * crossVector);
            double m_inv = 1/this->mass + 1/other.mass + (r1.cross(N_hat)).transpose() * this->inertia_tensor.inverse() * (r1.cross(N_hat)) + r2.cross(N_hat).squaredNorm() / I_CM;
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) + 2 * r2.cross(N_hat).dot(other.w_CM) - 2 * r1.cross(N_hat).dot(this->w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * (this->inertia_tensor.inverse() * r1.cross(N_hat));
            other.v_CM -= Impulse * N_hat / other.mass;
            other.w_CM -= Impulse * r2.cross(N_hat) / I_CM;
        }
        else {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r1 = this->atoms[i].coor - this->coor_CM;
            Eigen::Vector3d r2 = other.atoms[j].coor - other.coor_CM;
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(this->linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(this->linear_axis) * this->linear_axis).normalized();
            double I_CM_1 = crossVector.dot(this->inertia_tensor * crossVector);
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(other.linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(other.linear_axis) * other.linear_axis).normalized();
            double I_CM_2 = crossVector.dot(other.inertia_tensor * crossVector);
            double m_inv = 1/this->mass + 1/other.mass + (r1.cross(N_hat).squaredNorm() / I_CM_1) + (r2.cross(N_hat).squaredNorm() / I_CM_2);
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) + 2 * r2.cross(N_hat).dot(other.w_CM) - 2 * r1.cross(N_hat).dot(this->w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * r1.cross(N_hat) / I_CM_1;
            other.v_CM -= Impulse * N_hat / other.mass;
            other.w_CM -= Impulse * r2.cross(N_hat) / I_CM_2;
        }
    }
    else if (this->is_polyatomic) {
        if (!this->is_linear) {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r1 = this->atoms[i].coor - this->coor_CM;
            double m_inv = 1/this->mass + 1/other.mass + r1.cross(N_hat).transpose() * this->inertia_tensor.inverse() * (r1.cross(N_hat));
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) - 2 * r1.cross(N_hat).dot(this->w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * (this->inertia_tensor.inverse() * r1.cross(N_hat));
            other.v_CM -= Impulse * N_hat / other.mass;
        }
        else {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r1 = this->atoms[i].coor - this->coor_CM;
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(this->linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(this->linear_axis) * this->linear_axis).normalized();
            double I_CM = crossVector.dot(this->inertia_tensor * crossVector);
            double m_inv = 1/this->mass + 1/other.mass + (r1.cross(N_hat).squaredNorm() / I_CM);
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) - 2 * r1.cross(N_hat).dot(this->w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * r1.cross(N_hat) / I_CM;
            other.v_CM -= Impulse * N_hat / other.mass;
        }
    }
    else if (other.is_polyatomic) {
        if (!other.is_linear) {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r2 = other.atoms[j].coor - other.coor_CM;
            double m_inv = 1/this->mass + 1/other.mass + r2.cross(N_hat).transpose() * other.inertia_tensor.inverse() * (r2.cross(N_hat));
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) + 2 * r2.cross(N_hat).dot(other.w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            other.v_CM -= Impulse * N_hat / other.mass;
            other.w_CM -= Impulse * (other.inertia_tensor.inverse() * r2.cross(N_hat));
        }
        else {
            Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
            Eigen::Vector3d r2 = other.atoms[j].coor - other.coor_CM;
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(other.linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(other.linear_axis) * other.linear_axis).normalized();
            double I_CM = crossVector.dot(other.inertia_tensor * crossVector);
            double m_inv = 1/this->mass + 1/other.mass + (r2.cross(N_hat).squaredNorm() / I_CM);
            double vel = 2 * N_hat.dot(other.v_CM - this->v_CM) + 2 * r2.cross(N_hat).dot(other.w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            other.v_CM -= Impulse * N_hat / other.mass;
            other.w_CM -= Impulse * r2.cross(N_hat) / I_CM;
        }
    }
    else {
        Eigen::Vector3d N_hat = (this->atoms[i].coor - other.atoms[j].coor).normalized();
        double m_inv = 1/this->mass + 1/other.mass;
        double vel = 2 * N_hat.dot(other.v_CM - this->v_CM);
        double Impulse = vel / m_inv;
        this->v_CM += Impulse * N_hat / this->mass;
        other.v_CM -= Impulse * N_hat / other.mass;
    }
}
bool molecule::check_collision_wall() const {
    if (this->coor_CM(0) - this->radius <= 0 || this->coor_CM(0) + this->radius >= Constants::SCREEN(0) ||
        this->coor_CM(1) - this->radius <= 0 || this->coor_CM(1) + this->radius >= Constants::SCREEN(1) ||
        this->coor_CM(2) - this->radius <= 0 || this->coor_CM(2) + this->radius >= Constants::SCREEN(2)) {
        if (!this->is_polyatomic) {
            return true;
        }
        for (int i = 0; i < this->numatoms; i++) {
            if (this->atoms[i].checkCollisionWall() != 0) {
                return true;
            }
        }
    }
    return false;
}
void molecule::collision_wall() {
    if (this->coor_CM(0) - this->radius <= 0 || this->coor_CM(0) + this->radius >= Constants::SCREEN(0) ||
        this->coor_CM(1) - this->radius <= 0 || this->coor_CM(1) + this->radius >= Constants::SCREEN(1) ||
        this->coor_CM(2) - this->radius <= 0 || this->coor_CM(2) + this->radius >= Constants::SCREEN(2)) {
        for (int i = 0; i < this->numatoms; i++) {
            int a = this->atoms[i].checkCollisionWall();
            this->cal_collision_wall(i, a);
            if (a != 0) {
                break;
            }
        }
    }
}
void molecule::cal_collision_wall(const int& i, const int& num) {
    Eigen::Vector3d r = this->atoms[i].coor - this->coor_CM;
    Eigen::Vector3d N_hat;
    Eigen::Vector3d a;
    if (this->is_polyatomic) {
        if (!this->is_linear) {
            if (num == 0) {
                return;
            }
            else if (num == 1) {
                N_hat = Eigen::Vector3d(1.0,0.0,0.0);
                a = N_hat * (-this->atoms[i].coor(0) + this->atoms[i].radius);
            }
            else if (num == 2) {
                N_hat = Eigen::Vector3d(0.0,1.0,0.0);
                a = N_hat * (-this->atoms[i].coor(1) + this->atoms[i].radius);
            }
            else if (num == 3) {
                N_hat = Eigen::Vector3d(0.0,0.0,1.0);
                a = N_hat * (-this->atoms[i].coor(2) + this->atoms[i].radius);
            }
            else if (num == 4) {
                N_hat = Eigen::Vector3d(-1.0,0.0,0.0);
                a = N_hat * (this->atoms[i].coor(0) + this->atoms[i].radius - Constants::SCREEN(0));
            }
            else if (num == 5) {
                N_hat = Eigen::Vector3d(0.0,-1.0,0.0);
                a = N_hat * (this->atoms[i].coor(1) + this->atoms[i].radius - Constants::SCREEN(1));
            }
            else {
                N_hat = Eigen::Vector3d(0.0,0.0,-1.0);
                a = N_hat * (this->atoms[i].coor(2) + this->atoms[i].radius - Constants::SCREEN(2));
            }
            double m_inv = 1/this->mass + r.cross(N_hat).transpose() * this->inertia_tensor.inverse() * (r.cross(N_hat));
            double vel = -2 * N_hat.dot(this->v_CM) - 2 * r.cross(N_hat).dot(this->w_CM);
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * (this->inertia_tensor.inverse() * r.cross(N_hat));
        }
        else {
            if (num == 0) {
                return;
            }
            else if (num == 1) {
                N_hat = Eigen::Vector3d(1.0,0.0,0.0);
                a = N_hat * (-this->atoms[i].coor(0) + this->atoms[i].radius);
            }
            else if (num == 2) {
                N_hat = Eigen::Vector3d(0.0,1.0,0.0);
                a = N_hat * (-this->atoms[i].coor(1) + this->atoms[i].radius);
            }
            else if (num == 3) {
                N_hat = Eigen::Vector3d(0.0,0.0,1.0);
                a = N_hat * (-this->atoms[i].coor(2) + this->atoms[i].radius);
            }
            else if (num == 4) {
                N_hat = Eigen::Vector3d(-1.0,0.0,0.0);
                a = N_hat * (this->atoms[i].coor(0) + this->atoms[i].radius - Constants::SCREEN(0));
            }
            else if (num == 5) {
                N_hat = Eigen::Vector3d(0.0,-1.0,0.0);
                a = N_hat * (this->atoms[i].coor(1) + this->atoms[i].radius - Constants::SCREEN(1));
            }
            else {
                N_hat = Eigen::Vector3d(0.0,0.0,-1.0);
                a = N_hat * (this->atoms[i].coor(2) + this->atoms[i].radius - Constants::SCREEN(2));
            }
            Eigen::Vector3d crossVector;
            if (Eigen::Vector3d({1.0,0.0,0.0}).dot(this->linear_axis) > 0.99) {
                crossVector = Eigen::Vector3d(0.0, 1.0, 0.0);
            } else {
                crossVector = Eigen::Vector3d(1.0, 0.0, 0.0);
            }
            crossVector = (crossVector - crossVector.dot(this->linear_axis) * this->linear_axis).normalized();
            double I_CM = crossVector.dot(this->inertia_tensor * crossVector);
            double m_inv = 1/this->mass + r.cross(N_hat).squaredNorm() / I_CM;
            double vel = -2 * N_hat.dot(this->v_CM) - 2 * this->w_CM.dot(r.cross(N_hat));
            double Impulse = vel / m_inv;
            this->v_CM += Impulse * N_hat / this->mass;
            this->w_CM += Impulse * r.cross(N_hat) / I_CM;
        }
        this->coor_CM += a;
        for (int i = 0; i < this->numatoms; i++) {
            this->atoms[i].coor = this->atoms_coor[i] + this->coor_CM;
        }
    }
    else {
        if (num == 0) {
            return;
        }
        else if (num == 1) {
            this->v_CM(0) *= -1.0;
            a = Eigen::Vector3d(-this->atoms[i].coor(0) + this->atoms[i].radius, 0.0, 0.0);
        }
        else if (num == 2) {
            this->v_CM(1) *= -1.0;
            a = Eigen::Vector3d(0.0, -this->atoms[i].coor(1) + this->atoms[i].radius, 0.0);
        }
        else if (num == 3) {
            this->v_CM(2) *= -1.0;
            a = Eigen::Vector3d(0.0, 0.0, -this->atoms[i].coor(2) + this->atoms[i].radius);
        }
        else if (num == 4) {
            this->v_CM(0) *= -1.0;
            a = -Eigen::Vector3d(this->atoms[i].coor(0) + this->atoms[i].radius - Constants::SCREEN(0), 0.0, 0.0);
        }
        else if (num == 5) {
            this->v_CM(1) *= -1.0;
            a = -Eigen::Vector3d(0.0, this->atoms[i].coor(1) + this->atoms[i].radius - Constants::SCREEN(1), 0.0);
        }
        else {
            this->v_CM(2) *= -1.0;
            a = -Eigen::Vector3d(0.0, 0.0, this->atoms[i].coor(2) + this->atoms[i].radius - Constants::SCREEN(2));
        }
        this->coor_CM += a;
        this->atoms[0].coor = this->atoms_coor[0] + this->coor_CM;
    }
}
bool molecule::operator==(const molecule& rhs) const {
    if (this == &rhs) {
        return true;
    }
    else {
        return false;
    }
}
bool molecule::operator!=(const molecule& rhs) const {
    return !(*this == rhs);
}
void molecule::del(vector<molecule*>& Particles, vector<vector<vector<int>>>* checking_lst, int* num_molecules) {
    int idx_to_del = this->index;
    if (this->index == -1) {
        delete this;  // 현재 객체를 삭제
        return;
    }
    if (idx_to_del >= static_cast<int>(Particles.size()) || Particles[idx_to_del] != this) {
        cerr << "Error: Index mismatch in del()! Absorting deletion." << endl;
        return;
    }
    molecule* temp = Particles[Particles.size() - 1];
    Particles[Particles.size() - 1] = this;
    Particles[idx_to_del] = temp;
    this->index = temp->index;
    temp->index = idx_to_del;
    if (idx_to_del < static_cast<int>(Particles.size() - 1)) {
        for (size_t i = 0; i < Particles.size() - 1; i++) {
            if (static_cast<int>(i) < idx_to_del) {
                vector<int> temp_vec = (*checking_lst)[i][idx_to_del - i - 1];
                (*checking_lst)[i][idx_to_del - i - 1] = (*checking_lst)[i][Particles.size() - i - 2];
                (*checking_lst)[i][Particles.size() - i - 2] = temp_vec;
            }
            else if (static_cast<int>(i) > idx_to_del) {
                vector<int> temp_vec = (*checking_lst)[idx_to_del][i - idx_to_del - 1];
                (*checking_lst)[idx_to_del][i - idx_to_del - 1] = vector<int>{(*checking_lst)[i][Particles.size() - i - 2][1], (*checking_lst)[i][Particles.size() - i - 2][0]};
                (*checking_lst)[i][Particles.size() - i - 2] = vector<int>{temp_vec[1], temp_vec[0]};
            }
        }
    }
    for (int i = Particles.size() - 2; i >= 0; i--) {
        (*checking_lst)[i].pop_back();
    }
    (*checking_lst).pop_back(); 
    Particles.pop_back();
    *num_molecules -= 1;
    delete this;  // 현재 객체를 삭제
}
Eigen::Vector3d molecule::get_w_CM(const double& Energy, const Eigen::Vector3d& vector) const {
    if (!this->is_polyatomic) {
        return Eigen::Vector3d(0.0,0.0,0.0);
    }
    double I = vector.dot(this->inertia_tensor * vector);
    return sqrt(2 * Energy/ I) * vector;
}
Eigen::Vector3d molecule::get_v_CM(const double& Energy, const Eigen::Vector3d& vector) const {
    return sqrt(2 * Energy / this->mass) * vector;
}
string molecule::get_name() const {
    return this->name;
}
molecule::~molecule() {}
