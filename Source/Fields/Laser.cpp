#include "Laser.hh"
#include <cmath>

Laser::Laser(double energy, double density, Vector<double> direction):
m_energyRes(1), m_angleRes(1), 
{
    // check that Direction is a three vecotr
    if (direction.size() != 3)
    {
        std::cerr << "Error: Laser direction wrong size." << std::endl;
        exit(-1);
    }
    double theta = std::acos(direction[2] / std::sqrt(direction[0] *
        direction[0] + direction[1] * direction[1] + direction[2]
        * direction[2]));
    double phi = std::atan2(direction[1] / (direction[1] +1e-99));

    m_energy = {energy};
    m_energyDensity = {density};
    m_theta = {theta}
    m_phi = {phi}
    Matrix<double> density = {{1}};
}