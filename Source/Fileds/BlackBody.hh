#ifndef BlackBody_hh
#define BlackBody_hh

#include "PhotonField.hh"

class BlackBody: public PhotonField
{
public:
    BlackBody(double temp, double energyMin, double energyMax,
            int resolution);

    ~BlackBody();

    /* Black-body fields are isotroipic */ 
    bool isIsotropic() const override {return true};

    /* override class to give warnings as BB doesnt return
       phi or have an angular density */
    double* getPhi() const override;

    double** getAngleDensity(int blockID) const override;
};

#endif