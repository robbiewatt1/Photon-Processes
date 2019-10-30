#ifndef PhotonField_hh
#define PhotonField_hh


class PhotonField
{
public:
    /* defult constructor / destructor */
    PhotonField();

     virtual ~PhotonField() = 0;

    /* Following getter methods allow acces to the photon field 
       energy, angle and density variables. */
    virtual int getEnergyRes() const {return m_energyRes;}

    virtual int getAngleRes() const {return m_angleRes;}

    virtual double* getEnergy() const {return m_energyAxis;}

    virtual double* getEnergyDensity() const {return m_energySpec;}

    virtual double* getTheta() const {return m_thetaAxis;}

    virtual double* getPhi() const {return m_phiAxis;}

    virtual double** getAngleDensity(int blockID) const 
            {return m_anguleDensity[blockID];}

    virtual double getBlockDens(int blockID) const {return
        m_blockDensity[blockID];}

    virtual int getBlockID(int pos[3]) const = 0;

    virtual void getBlockDims(int* dims) const = 0;

protected:

    double* m_energyAxis;
    double* m_energySpec;
    double* m_thetaAxis;
    double* m_phiAxis;

    double** m_anguleDensity;
    int* m_blockPos;
    double m_blockDensity;

    int m_nBlocks;
    int m_energyRes;
    int m_angleRes;
};

#endif