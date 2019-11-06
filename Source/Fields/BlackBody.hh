#ifndef BlackBody_hh
#define BlackBody_hh

#include "PhotonField.hh"

class BlackBody: public PhotonField
{
public:
    BlackBody(double temp, double energyMin, double energyMax,
            int energyRes, int angularRes);

    ~BlackBody();

    /* Black-body fields are isotroipic */ 
    bool isIsotropic() const override {return true;}

    int fieldDimensions() const override {return 2;}

    /* override to give warnings as BB doesnt return
       phi or have an angular density */
    const Vector<double>& getPhi() const override;

    const Matrix<double>& getAngleDensity(int blockID) const override;
};

#endif