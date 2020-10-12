#include "PhotonScatter.hh"
#include "Numerics.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include <sys/stat.h>
#include "G4AutoLock.hh"

namespace { G4Mutex hdf5Mutex = G4MUTEX_INITIALIZER; }


PhotonScatter::PhotonScatter(PhotonField* field, double comMin):
PhotonProcess(field, comMin, "PhotonScatter")
{
    loadCrossSection();
}

#ifdef USEGP
PhotonScatter::PhotonScatter(PhotonField* field, int trainSize,
    double errorMax, double comMin, std::string saveDir):
PhotonProcess(field, comMin, trainSize, errorMax, saveDir, "PhotonScatter")
{
    loadCrossSection();
}

PhotonScatter::PhotonScatter(PhotonField* field, const G4String& gpDir,
    double errorMax, double comMin, std::string saveDir):
PhotonProcess(field, comMin, gpDir, errorMax,
    saveDir, "PhotonScatter")
{
    loadCrossSection();
}
#endif

PhotonScatter::~PhotonScatter()
{
}

G4VParticleChange* PhotonScatter::PostStepDoIt(const G4Track& aTrack,
        const G4Step& aStep)
{
   /* Get the properties of both the interacting gamma and a sampled 
       photon from the photon field */
    aParticleChange.Initialize(aTrack);
    const G4DynamicParticle *aDynamicGamma = aTrack.GetDynamicParticle();
    G4Material* aMaterial = aTrack.GetMaterial();
    int blockID = aMaterial->GetMaterialPropertiesTable()
        ->GetConstProperty("BlockID");

    /* Gamma properties */
    double gammaEnergy = aDynamicGamma->GetKineticEnergy();
    G4ThreeVector gammaDirection = aDynamicGamma->GetMomentumDirection();

    /* Rotations to gamma frame */
    G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1).cross(-gammaDirection);
    double rotationAngle = gammaDirection.angle(G4ThreeVector(0, 0, 1));
    m_rotaion  = G4RotationMatrix(rotationAxis, rotationAngle);
    G4RotationMatrix rotateBack = G4RotationMatrix(rotationAxis,
        -rotationAngle);

    /* Photon properties */
    double photonEnergy, comEnergy, photonTheta, photonPhi;
    samplePhotonField(blockID, gammaEnergy, photonEnergy, comEnergy, photonPhi);
    photonTheta = centreOfMassTheta(comEnergy, gammaEnergy, photonEnergy);

    /* Find particles properties in the COM frame*/
    double pairEnergy = electron_mass_c2 * std::sqrt(comEnergy) / 2.0;
    double pairPhi = 2.0 * pi * G4UniformRand();
    double pairTheta = samplePairAngle(comEnergy);

    G4LorentzVector gamma1Vector = G4LorentzVector(
        std::sin(pairTheta) * std::cos(pairPhi) * pairEnergy,
        std::sin(pairTheta) * std::sin(pairPhi) * pairEnergy,
        std::cos(pairTheta) * pairEnergy, pairEnergy);
    G4LorentzVector gamma2Vector = G4LorentzVector(
        -std::sin(pairTheta) * std::cos(pairPhi) * pairEnergy,
        -std::sin(pairTheta) * std::sin(pairPhi) * pairEnergy,
        -std::cos(pairTheta) * pairEnergy, pairEnergy);

    double beta = std::sqrt((gammaEnergy * gammaEnergy + photonEnergy
        * photonEnergy + 2.0 * gammaEnergy * photonEnergy
        * std::cos(photonTheta)) / (gammaEnergy * gammaEnergy + photonEnergy
            * photonEnergy + 2.0 * gammaEnergy * photonEnergy));
    gamma1Vector.boostZ(beta);
    gamma2Vector.boostZ(beta);

    /* Apply rotation into frame with gamma along z axis */  
    double thetaGamma = std::atan2(photonEnergy * std::sin(photonTheta),
                 gammaEnergy + photonEnergy * std::cos(photonTheta));
    gamma1Vector.rotateY(thetaGamma);
    gamma2Vector.rotateY(thetaGamma);

    /* Apply rotation around the gamma axis by -photonPhi */
    gamma1Vector.rotateZ(-photonPhi);
    gamma2Vector.rotateZ(-photonPhi);

    /* Apply final rotation into simulation frame */
    gamma1Vector = rotateBack(gamma1Vector);
    gamma2Vector = rotateBack(gamma2Vector);

    /* Add new electron and positron and kill the gamma ray */
    G4DynamicParticle* photon1 = new G4DynamicParticle(
            G4Gamma::Gamma(), gamma1Vector);            
    G4DynamicParticle* photon2 = new G4DynamicParticle(
            G4Gamma::Gamma(), gamma2Vector);
    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(photon1);
    aParticleChange.AddSecondary(photon2);
    aParticleChange.ProposeMomentumDirection(G4ThreeVector(0,0,0));
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4bool PhotonScatter::IsApplicable(const G4ParticleDefinition&
        particle) const
{
    return (&particle == G4Gamma::Gamma());
}

double PhotonScatter::crossSection(double comEnergy) const
{
    return Numerics::interpolate1D(m_comEnergy, m_totalCrossSection,
        comEnergy);
}

double PhotonScatter::centreOfMassEnergy(double dynamicEnergy,
    double staticEnergy, double theta) const
{
    return 2.0 * dynamicEnergy * staticEnergy * (1.0 - std::cos(theta))
        / (electron_mass_c2 * electron_mass_c2);
}

double PhotonScatter::centreOfMassStatic(double comEnergy,
    double dynamicEnergy, double theta) const
{
    return electron_mass_c2 * electron_mass_c2 * comEnergy
        / (2.0 * dynamicEnergy * (1.0 - std::cos(theta)));
}

double PhotonScatter::centreOfMassDynamic(double comEnergy,
        double staticEnergy, double theta) const
{
    return electron_mass_c2 * electron_mass_c2 * comEnergy
        / (2.0 * staticEnergy * (1.0 - std::cos(theta)));
}

double PhotonScatter::centreOfMassTheta(double comEnergy,
    double dynamicEnergy, double staticEnergy) const
{
    return std::acos(1.0 - electron_mass_c2 * electron_mass_c2
        * comEnergy / (2.0 * dynamicEnergy * staticEnergy));
}

void PhotonScatter::loadCrossSection()
{
    G4AutoLock lock(&hdf5Mutex);

    std::string fileName = theProcessName + "_total.h5";
    // Check if new file exisits in current dir
    struct stat buffer;
    if (stat (fileName.c_str(), &buffer) == 0)
    {
        // Load file from run dir
        m_comEnergy.open(fileName, "/s");
        m_totalCrossSection.open(fileName, "/sigma");
    } else
    {
        // Load file from Install
        try
        {
            std::string dataDir(getenv("Photon_Process_Data_Dir"));
            std::string filePath = dataDir + "/" + fileName;
            m_comEnergy.open(filePath, "/s");
            m_totalCrossSection.open(filePath, "/sigma");
        } catch (const std::exception& e)
        {
            std::cout << e.what() << std::endl;
            std::cout << "You probably haven't run PhotonProcess.sh"
                << std::endl;
            std::abort();
        }

    }
}