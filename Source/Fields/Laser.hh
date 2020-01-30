#ifndef Laser_hh
#define Laser_hh

#include "PhotonField.hh"
#include "Vector.hh"


class Laser: public PhotonField
{
public:
    /* Defult constructor for Laser radiation field.
    energy: Photon energy in MeV
    density: Photon density in mm^-3
    direction: Direction of laser
    */
    explicit Laser(double energy, double density, Vector<double> direction);

    /* Laser fields are not isotroipic */ 
    bool isIsotropic() const override {return false;}

    /* returns the dimensionality of the field. */
    int fieldDimensions() const override {return 2;}
};
#endif