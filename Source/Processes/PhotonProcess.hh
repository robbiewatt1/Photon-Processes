#ifndef PhotonProcess_hh
#define PhotonProcess_hh

#include "G4VDiscreteProcess.hh"
#include "PhotonField.hh"

#ifdef USEGP
    #include "GaussianProcess.hh"
#endif


class PhotonProcess: public G4VDiscreteProcess
{
public:
    explicit PhotonProcess(PhotonField* field, double comMin,
        const G4String& name, G4ProcessType type = fUserDefined);

#ifdef USEGP
    explicit PhotonProcess(PhotonField* field, double comMin,
        int trainSize, double errorMax, bool save,
        const G4String& name, G4ProcessType type = fUserDefined);

    explicit PhotonProcess(PhotonField* field,
        double comMin, const G4String& gpDir,
        int trainSize, double errorMax, bool save,
        const G4String& name,
        G4ProcessType type = fUserDefined);
#endif

    virtual ~PhotonProcess();

    /* Method to set the paramters used by the gaussian process */
    void setParamsGP(bool save, int trainSize, double errorMax);

    /* Method to load a pre-saved GP from a file */
    void loadGP(std::string gpDir);

    /* Method to caculate the mean free path for interacting gamma. */
    G4double GetMeanFreePath(const G4Track& track, G4double, 
            G4ForceCondition*) override;

    /* Method to add new particles and edit interacting ones */
    virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
            const G4Step& aStep) = 0;

/* The following methods must be overriden for a new process */
protected:

    /* Returns the total cross-section for the process */
    virtual double crossSection(double comEnergy) const = 0;

    /* Returns the differential cross section of the process */
    virtual double diffCrossSection(double comEnergy, double theta) const = 0;

    /* Calculates the centre of mass energy. dynamicEnergy is the 
       interacting particle and static energy is a photon from the field */
    virtual double centreOfMass(double dynamicEnergy, double staticEnergy,
        double theta) const = 0;

protected:
    /* Samples the energy, angle from the photon field */
    void samplePhotonField(int blockID, double gammaEnergy, double& photonEnergy,
        double& comEnergy, double& photonPhi);

    /* Sampeles a scattering angle from the differential cross-section.
       Assumes scattering is maximum on axis.  */
    double samplePairAngle(double comEnergy);

    /* returns the values of theta and phi in the rotated frame.
       angleIn = [theta, phi] angleOut = [theta, phi]*/
    void rotateThetaPhi(double angleIn[2], double angleOut[2],
        const G4RotationMatrix& rotaion);

private:
    PhotonField* m_field;      // static X ray field
    double m_comMin;           // Min COM energy of interest
    double m_minDensity;       // Min angula density of interest
    G4RotationMatrix m_rotaion; // Rotation matrix to gamma frame

#ifdef USEGP
    double m_errorMax;       // Max eroor before using GP
    bool m_save;             // Bool deciding if GP is saved
    GaussianProcess* m_gp;   // GP class handling the Gaussian Process
#endif

};

#endif