#ifndef PhotonProcess_hh
#define PhotonProcess_hh

#include "G4VDiscreteProcess.hh"
#include "PhotonField.hh"


class PhotonProcess: public G4VDiscreteProcess
{
public:
    explicit PhotonProcess(PhotonField* field,
        int trainSize, double errorMax, bool save,
        const G4String& name, G4ProcessType type = fUserDefined);

    explicit BreitWheelerGP(PhotonField* field,
        const G4String& gpDir,
        int trainSize, double errorMax, bool save,
        const G4String& name = "BreitWheelerGP",
        G4ProcessType type = fUserDefined);

    ~BreitWheelerGP();

    /* Method to set the paramters used by the gaussian process */
    void setParamsGP(bool save, int trainSize, double errorMax);

    /* Method to load a pre-saved GP from a file */
    void loadGP(std::string gpDir;);

    /* Method to caculate the mean free path for interacting gamma. */
    G4double GetMeanFreePath(const G4Track& track, G4double, 
            G4ForceCondition*) override;

    virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
            const G4Step& aStep) = 0;

/* The following methods must be overriden for a new process */
protected:

    /* Returns the total cross-section for the process */
    double crossSection(double comEnergy) const = 0;

    /* Returns the differential cross section of the process */
    double diffCrossSection(double comEnergy, double theta) const = 0;

protected:
    /* Samples the energy, angle from the photon field */
    void SamplePhotonField(PhotonField* field, double gammaEnergy,
            double& photonEnergy, double& comEnergy, double& photonPhi);
};

#endif