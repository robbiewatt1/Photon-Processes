#include "Laser.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

Laser::Laser(double energy, double energyDensity, Vector<double> direction)
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
    double phi = std::atan2(direction[1], (direction[0] + 1e-99));

    m_nBlocks = 1;
    m_energyRes = 1;
    m_angleRes = 1;
    m_energy = {energy};
    m_energyDensity = {energyDensity};
    m_theta = {theta};
    m_phi = {phi};
    m_angleDensity = {{{1.0}}};
}