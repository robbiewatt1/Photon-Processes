#ifndef PhotonScatter_hh
#define PhotonScatter_hh

#include "PhotonProcess.hh"

class PhotonScatter: public PhotonProcess
{
public:
    /* Defult no GP constructor */
    explicit PhotonScatter(PhotonField* field, const std::string& dataFile,
        double comMin);

#ifdef USEGP
    /* GP constructor new */
    explicit PhotonScatter(PhotonField* field, const std::string& dataFile,
        int trainSize, double errorMax, double comMin,
        std::string saveDir = "");

    /* GP constructor using previously trained GP */
    explicit PhotonScatter(PhotonField* field, const std::string& dataFile,
        const G4String& gpDir, double errorMax, double comMin,
        std::string saveDir = "");
#endif

    ~PhotonScatter();  

    G4bool IsApplicable(const G4ParticleDefinition&
        particle) const override;
    
    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
        const G4Step& aStep) override;

public:
    /* Returns the total cross-section for the process */
    double crossSection(double comEnergy) const override;

    /* Calculates the centre of mass energy. dynamicEnergy is the 
       interacting particle and static energy is a photon from the field */
    double centreOfMassEnergy(double dynamicEnergy,
        double staticEnergy, double theta) const override;

    /* Same equation as centreOfMass but rearanged so static energy is
       returned */
    double centreOfMassStatic(double comEnergy,
        double dynamicEnergy, double theta) const override;

    /* Same equation as centreOfMass but rearanged so dynamic energy is
       returned */
    double centreOfMassDynamic(double comEnergy,
        double staticEnergy, double theta) const override;

    /* Same equation as centreOfMass but rearanged so angle between particles
       is returned */
    double centreOfMassTheta(double comEnergy, double dynamicEnergy,
        double staticEnergy) const override;

private:
    /* Opens the file containing the differential
       cross-section data */
    void openDataFile(const std::string& fileName);

    Vector<double> m_comEnergy;
    Vector<double> m_comTheta;
    Vector<double> m_totalCrossSection;
};
#endif