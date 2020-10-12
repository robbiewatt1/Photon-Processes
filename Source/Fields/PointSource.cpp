#include "PointSource.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>
#include "G4Threading.hh"
#include <stdlib.h>

#include "G4AutoLock.hh"
namespace { G4Mutex hdf5Mutex = G4MUTEX_INITIALIZER; }


PointSource::PointSource(const Vector<double>& res,
    const Vector<double>& extent)
{
    // check that Direction is a three vecotr
    if (res.size() != 3 || extent.size() != 3)
    {
        std::cerr << "Error: Point source extent and resolution." << std::endl;
        exit(-1);
    }

    m_nBlocks = res[0] * res[1] * res[2];
    m_sourceDistance2 = Vector<double>(m_nBlocks);
    m_sourceTheta = Vector<double>(m_nBlocks);
    m_sourcePhi = Vector<double>(m_nBlocks);

    int index = 0;
    for (int i = 0; i < res[0]; ++i)
    {
        for (int j = 0; j < res[1]; ++j)
        {
            for (int k = 0; k < res[2]; ++k)
            {
                double x = (0.5 + i - res[0] / 2.0) * extent[0] / res[0];
                double y = (0.5 + j) * extent[1] / res[1];
                double z = (0.5 + k - res[2] / 2.0) * extent[2] / res[2];
                m_sourceDistance2[index] = x * x + y * y + z * z;
                m_sourcePhi[index] = std::atan(y / x);
                m_sourceTheta[index] = std::acos(z
                    / std::sqrt(m_sourceDistance2[index]));    
                index++;   
            }
        }
    }
    m_angleRes = 1;
}

PointSource::~PointSource()
{
}

void PointSource::setSpectrumGaussian(double meanEnergy, double sigEnergy,
    double density, int energyRes, double energyMin, double energyMax)
{
    m_energyRes = energyRes;
    m_energy  = Vector<double>(m_energyRes);
    m_energyDensity = Vector<double>(m_energyRes);

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

    m_energyDensity = Vector<double>(m_energyRes);
    for (int i = 0; i < m_energyRes; ++i)
    {
        m_energyDensity[i] = density / (sigEnergy * std::sqrt(pi))
            * std::exp(-(m_energy[i] - meanEnergy) * (m_energy[i]
                - meanEnergy) / (sigEnergy * sigEnergy));
    }
}

void PointSource::setSpectrumFile(std::string fileName)
{
    G4AutoLock lock(&hdf5Mutex);

    m_file = new H5::H5File(fileName, H5F_ACC_RDWR);
    m_energy.open(fileName, "/Energy/Axis");
    m_energyDensity.open(fileName, "/Energy/Spectrum");
    m_energyRes = m_energy.size();
}

const Matrix<double>& PointSource::getAngleDensity(int blockID)
{
    m_angleDensity = {{{1.0 / (m_sourceDistance2[blockID])}}};
    return m_angleDensity[0];
}

const Vector<double>& PointSource::getTheta(int blockID)
{
    m_theta = {m_sourceTheta[blockID]};
    return m_theta;
}

const Vector<double>& PointSource::getPhi(int blockID)
{
    m_phi = {m_sourcePhi[blockID]};
    return m_phi;
}