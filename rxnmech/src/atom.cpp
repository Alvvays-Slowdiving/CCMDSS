#include "atom.hpp"

atom::atom() : radius(0), mass(0), coor(Eigen::Vector3d::Zero()), name("") {}

atom::atom(double r, double m, const Eigen::Vector3d& coord, const std::string& n)
    : radius(r), mass(m), coor(coord), name(n) {}

bool atom::checkCollision(const atom& other) const {
    return (coor - other.coor).norm() <= (radius + other.radius);
}

int atom::checkCollisionWall(const Eigen::Vector3d& screen) const {
    if (coor[0] - radius <= 0)
        return 1;
    else if (coor[1] - radius <= 0)
        return 2;
    else if (coor[2] - radius <= 0)
        return 3;
    else if (coor[0] >= screen[0] - radius)
        return 4;
    else if (coor[1] >= screen[1] - radius)
        return 5;
    else if (coor[2] >= screen[2] - radius)
        return 6;
    else
        return 0;
}
std::string atom::toString() const {
    return name;
}