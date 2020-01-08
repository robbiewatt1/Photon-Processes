#ifndef BlackBody_hh
#define BlackBody_hh

#include "PhotonField.hh"

class BlackBody: public PhotonField
{
public:
    /* Defult constructor for Blackbody radiation field.
    temp: Temperature of field in MeV.
    energyMin: Lowest energy considered in MeV.
    energyMax: Maximum energy considered in MeV
    energyRes: Resolusion for energy integration.
    angularRes: Resolusion for angular integration.
    */
    explicit BlackBody(double temp, double energyMin, double energyMax,
            int energyRes, int angularRes);

    
    ~BlackBody();

    /* Black-body fields are isotroipic */ 
    bool isIsotropic() const override {return true;}

    /* returns the dimensionality of the field. */
    int fieldDimensions() const override {return 1;}

    /* Override to give warnings as BB doesnt return
       phi or have an angular density */
    const Vector<double>& getPhi() const override;

    /* Override to give warnings as BB doesnt return
       an angular density */
    const Matrix<double>& getAngleDensity(int blockID) const override;
};

#endif