#include "Gaussian.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

Gaussian::Gaussian(double meanEnergy, double sigEnergy, double density,
        int energyRes, Vector<double> direction, double energyMin,
        double energyMax)
{
    // check that Direction is a three vecotr
    if (direction.size() != 3)
    {
        std::cerr << "Error: Laser direction wrong size." << std::endl;
        exit(-1);
    }
    m_nBlocks = 1;
    m_energyRes = energyRes;
    m_energy  = Vector<double>(m_energyRes);


    // If no min / max provided, use 4 s. d. from mean
    if (energyMin < 0)
    {
        energyMin = meanEnergy - 4.0 * sigEnergy;
        energyMin = energyMin < 0 ? 0 : energyMin;
    }
    if (energyMax < 0)
    {
        energyMax = meanEnergy + 4.0 * sigEnergy;
    }

    double energyDelta = (energyMax - energyMin) / m_energyRes;
    for (int i = 0; i < m_energyRes; i++) {m_energy[i] = i * energyDelta
            + energyMin;}

    double theta = std::acos(direction[2] / std::sqrt(direction[0] *
        direction[0] + direction[1] * direction[1] + direction[2]
        * direction[2]));
    double phi = std::atan2(direction[1], (direction[1] + 1e-99));

    m_energyDensity = Vector<double>(m_energyRes);
    for (int i = 0; i < m_energyRes; ++i)
    {
        m_energyDensity[i] = 1.0 / (sigEnergy * std::sqrt(pi))
            * std::exp(-(m_energy[i] - meanEnergy) * (m_energy[i]
                - meanEnergy) / (sigEnergy * sigEnergy));
    }

    m_nBlocks = 1;
    m_angleRes = 1;
    m_theta = {theta};
    m_phi = {phi};
    m_angleDensity = {{{1.0}}};
}

