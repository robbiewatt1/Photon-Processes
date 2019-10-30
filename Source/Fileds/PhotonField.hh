#ifndef PhotonField_hh
#define PhotonField_hh


class PhotonField
{
public:

    /* Override by returning true for isotropic or false 
       for nonisotroipic */
    virtual bool isIsotropic() const = 0;

    /* Following getter methods allow acces to the photon field 
       energy, angle and density variables. */
    virtual int getEnergyRes() const {return m_energyRes;}

    virtual int getAngleRes() const {return m_angleRes;}

    virtual double* getEnergy() const {return m_energy;}

    virtual double* getEnergyDensity() const {return m_energyDensity;}

    virtual double* getTheta() const {return m_theta;}

    virtual double* getPhi() const {return m_phi;}

    virtual double** getAngleDensity() const {return m_angleDensity;}

    virtual int getBlockID(int pos[3]) const {return 0;}


protected:

    int m_energyRes;
    int m_angleRes;
    double* m_energy;
    double* m_energyDensity;
    double* m_theta;
    double* m_phi;
    double** m_angleDensity;
//    int* m_blockPos;
//    double m_blockDensity;

};

#endif