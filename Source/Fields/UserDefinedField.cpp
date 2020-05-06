#include "UserDefinedField.hh"
#include "H5Cpp.h"

UserDefinedField::UserDefinedField(const std::string fileName):
m_fileName(fileName)
{
    // Check number of blocks in file
    file = new H5::H5File(fileName, H5F_ACC_RDWR);
    
    m_energy.open(fileName, "/EnergyAxis");
    m_energyDensity.open(fileName, "/EnergyDensity");
    m_theta.open(fileName, "/ThetaAxis");
    m_phi.open(fileName, "/PhiAxis");
    m_energyRes = m_energy.size();

    // Check that theta and pho are same size
    if(m_phi.size() != m_theta.size())
    {
        std::cerr << "Error: Different resolutions for theta and phi is not"
                  << "currently supported." << std::endl;
        std::exit(-1);
    }
    m_angleRes = m_theta.size();


    H5::Group* group = new H5::Group(file->openGroup("/AngleDensity"));
    hsize_t blocks;
    H5Gget_num_objs(group->getId(), &blocks);
    m_nBlocks = int(blocks);
    m_blockLocations = Vector<Vector<double>>(m_nBlocks);
    m_blocksize  = Vector<Vector<double>>(m_nBlocks);
    m_angleDensity = Vector<Matrix<double>>(m_nBlocks);

    // Loop over volumes
    for(int i = 0; i < m_nBlocks; i++)
    {
        std::string groupName = "/AngleDensity/" + std::to_string(i);
        m_angleDensity[i].open(fileName, groupName + "/density");
        m_blockLocations[i].open(fileName, groupName + "/location");
        m_blocksize[i].open(fileName, groupName + "/size");
    }
    delete group;
}

UserDefinedField::~UserDefinedField()
{
    delete file;
}

std::pair<Vector<Vector<double>>, Vector<Vector<double>>>
UserDefinedField::getDimensions() const
{
    return std::pair<Vector<Vector<double>>,
        Vector<Vector<double>>>(m_blockLocations, m_blocksize);
}
