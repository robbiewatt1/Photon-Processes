#ifndef PhotonField_hh
#define PhotonField_hh

#include "Matrix.hh"
#include "Vector.hh"

class PhotonField
{
public:

    /* Override by returning true for isotropic or false 
       for nonisotroipic */
    virtual bool isIsotropic() const = 0;

    /* returns the dimensionality of the field. i.e 2 for isotropic
       or 3 for nonisotropic */
    virtual int fieldDimensions() const = 0;

    /* Following getter methods allow acces to the photon field 
       energy, angle and density variables. */
    virtual int getEnergyRes() const {return m_energyRes;}

    virtual int getAngleRes() const {return m_angleRes;}

    virtual const Vector<double>& getEnergy() const {return m_energy;}

    virtual const Vector<double>& getEnergyDensity() const {
        return m_energyDensity;}

    virtual const Vector<double>& getTheta() const {return m_theta;}

    virtual const Vector<double>& getPhi() const {return m_phi;}

    virtual const Matrix<double>& getAngleDensity(int blockID) const
        {return m_angleDensity[blockID];}

    virtual int getBlockID(int pos[3]) const {return 0;}

    virtual int getNumBlocks() const {return 1.0;} 

protected:

    int m_energyRes;
    int m_angleRes;
    Vector<double> m_energy;
    Vector<double> m_energyDensity;
    Vector<double> m_theta;
    Vector<double> m_phi;
    Vector<Matrix<double>> m_angleDensity;
};

#endif