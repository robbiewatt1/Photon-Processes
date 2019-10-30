#include "BlackBody.hh"


BlackBody::BlackBody(double temp, double energyMin, double energyMax,
    int energyRes, int angularRes):
m_energyRes(energyRes), m_angleRes(angularRes)
m_anguleDensity(nullptr)
{
    m_energy  = new double [m_energyRes];
    m_theta   = new double [m_angleRes];

    double energyDelta = (energyMax - energyMin) / m_energyRes;
    double thetaDelta  = 2.0 * pi / m_angleRes;
    for (int i = 0; i < m_resolution; i++)
    {
        m_theta[i]   = i * thetaDelta;
        m_energy[i]  = i * energyDelta + energyMin;
    }

    m_energyDensity = new double [m_energyRes];
    for (int i = 0; i < count; ++i)
    {
        m_energyDensity[i] = 1.0 / (hbar_Planck * hbar_Planck * hbar_Planck
            * c_light * c_light * c_light * pi * pi) * m_energy[i]
            * m_energy[i] / std::exp(m_energy[i] / temp - 1.0);
    }
}

BlackBody::~BlackBody()
{
    delete [] m_energy;
    delete [] m_theta;
    delete [] m_energyDensity;
}
double* BlackBody::getPhi() const
{
    std::cerr << "Error: BlackBody is isotropic and should not call getPhi()."
        << std::endl;
    exit(1);
}

double** getAngleDensity(int blockID) const
{
    std::cerr << "Error: BlackBody is isotropic and should not call
        getAngleDensity()." << std::endl;
    exit(1);
}