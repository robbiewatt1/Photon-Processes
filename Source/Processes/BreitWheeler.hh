#ifndef BreitWheeler_hh
#define BreitWheeler_hh

#include "PhotonProcess.hh"

class BreitWheeler: public PhotonProcess
{
public:
    explicit BreitWheeler(PhotonField* field, double comMin = 4.0);

#ifdef USEGP
    explicit BreitWheeler(PhotonField* field, int trainSize, double errorMax,
        double comMin = 4.0, std::string saveDir = "");

    explicit BreitWheeler(PhotonField* field, const G4String& gpDir,
        double errorMax, double comMin = 4.0, std::string saveDir = "");
#endif

    ~BreitWheeler();

    G4bool IsApplicable(const G4ParticleDefinition& particle) const override;
    
    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
        const G4Step& aStep) override;

protected:
    /* Returns the total cross-section for the process */
    double crossSection(double comEnergy) const override;

    /* Returns the differential cross section of the process */
    double diffCrossSection(double comEnergy, double theta) const override;

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
};
#endif