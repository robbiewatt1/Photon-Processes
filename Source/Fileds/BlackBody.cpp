#include "BlackBody.hh"


BlackBody::BlackBody(double temp, double energyMin, double energyMax,
    int resolution)
{
    m_energy  = new double [m_resolution];
    m_theta   = new double [m_resolution];
}

BlackBody::BlackBody()
{
}