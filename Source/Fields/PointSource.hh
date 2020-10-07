#ifndef PointSource_hh
#define PointSource_hh

#include "PhotonField.hh"
#include <string>

class PointSource: public PhotonField
{
public:
    /* Constructor of class.file should have energy / density data in 
    MeV / 1 / (MeV mm^3) taken at 1mm above point source */
    PointSource(Vector<double> res, Vector<double> extent);
    
    ~PointSource();

    void setSpectrumGaussian(double meanEnergy, double sigEnergy,
        double density, int energyRes, double energyMin = -1,
        double energyMax = -1);

    void setSpectrumFile(std::string fileName);

    /* Black-body fields are isotroipic */ 
    bool isIsotropic() const override {return false;}

    /* returns the dimensionality of the field. */
    int fieldDimensions() const override {return 3;}
    
    /* Override to return density from point source (1/r^2) */
    const Matrix<double>& getAngleDensity(int blockID) override;

    const Vector<double>& getTheta(int blockID) override;

    const Vector<double>& getPhi(int blockID) override;


public:
    H5::H5File* m_file;
    Vector<double> m_sourceDistance2;
    Vector<double> m_sourceTheta;
    Vector<double> m_sourcePhi;

};
#endif