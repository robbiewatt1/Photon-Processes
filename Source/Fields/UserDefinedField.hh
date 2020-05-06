#ifndef UserDefinedField_hh
#define UserDefinedField_hh

#include "G4PVPlacement.hh"
#include "PhotonField.hh"
#include <utility>
  

class UserDefinedField: public PhotonField
{
public:
    /* Constructor of class. Must pass file name of 
       hdf5 file to be loaded including path */
    UserDefinedField(const std::string fileName);

    ~UserDefinedField();
    
    /* Override by returning true for isotropic or false 
       for nonisotroipic */
    bool isIsotropic() const override {return false;}

    /* returns the dimensionality of the field. i.e 2 for isotropic
       or 3 for nonisotropic */
    int fieldDimensions() const override {return 3;}

    /* Checks the photon field array is consistent with the detector setup */
    std::pair<Vector<Vector<double>>, Vector<Vector<double>>> getDimensions()
        const;

private:
    std::string m_fileName;
    H5::H5File* file;
    Vector<Vector<double>> m_blockLocations;
    Vector<Vector<double>> m_blocksize;
};
#endif