#ifndef Gaussian_hh
#define Gaussian_hh

#include "PhotonField.hh"

class Gaussian: public PhotonField
{
public:
    /* Defult constructor for Gaussian radiation field.
    meanEnergy: Mean photon energy in MeV
    sigEnergy: Variance of photon energy in MeV
    density: Photon density in mm^-3
    energyRes: Resolusion for energy integration.
    direction: Direction of laser
    */
    explicit Gaussian(double meanEnergy, double sigEnergy, double density,
        int energyRes, Vector<double> direction, double energyMin = -1,
        double energyMax = -1);
    
    ~Gaussian();

    /* Black-body fields are isotroipic */ 
    bool isIsotropic() const override {return true;}

    /* returns the dimensionality of the field. */
    int fieldDimensions() const override {return 3;}
};
#endif