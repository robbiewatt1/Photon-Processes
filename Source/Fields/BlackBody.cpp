#include "BlackBody.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <cmath>

BlackBody::BlackBody(double temp, double energyMin, double energyMax,
    int energyRes, int angularRes)
{
    m_nBlocks = 1;
    m_energyRes = energyRes;
    m_angleRes = angularRes;
    m_energy  = Vector<double>(m_energyRes);
    m_theta   = Vector<double>(m_angleRes);

    double energyDelta = (energyMax - energyMin) / m_energyRes;
    double thetaDelta  = CLHEP::pi / m_angleRes;
    for (int i = 0; i < m_angleRes; i++) m_theta[i] = i * thetaDelta;
    for (int i = 0; i < m_energyRes; i++) {m_energy[i] = i * energyDelta
            + energyMin;}

    m_energyDensity = Vector<double>(m_energyRes);
    for (int i = 0; i < m_energyRes; ++i)
    {
        m_energyDensity[i] = 1.0 / (hbar_Planck * hbar_Planck * hbar_Planck
            * c_light * c_light * c_light * pi * pi) * m_energy[i]
            * m_energy[i] / (std::exp(m_energy[i] / temp) - 1.0);
    }
}

BlackBody::~BlackBody()
{
}

const Vector<double>& BlackBody::getPhi(int blockID)
{
    std::cerr << "Error: BlackBody is isotropic and should not call getPhi()."
        << std::endl;
    exit(1);
}

const Matrix<double>& BlackBody::getAngleDensity(int blockID)
{
    std::cerr << "Error: BlackBody is isotropic and should not call"
        "getAngleDensity()." << std::endl;
    exit(1);
}