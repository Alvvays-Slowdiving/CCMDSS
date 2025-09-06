#ifndef ATOM_HPP
#define ATOM_HPP

#include "constants.hpp" // For Eigen::Vector3d and other constants

class atom {
public:
    double radius;
    double mass;
    Eigen::Vector3d coor;
    std::string name;
    atom(); // Default constructor
    atom(double r, double m, const Eigen::Vector3d& coord, const std::string& n);

    bool checkCollision(const atom& other) const;
    int checkCollisionWall(const Eigen::Vector3d& screen = Constants::SCREEN) const;
    std::string toString() const;
    ~atom() = default; // Default destructor
};

#endif // ATOM_HPP
