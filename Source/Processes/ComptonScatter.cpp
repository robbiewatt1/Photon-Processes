#include "ComptonScatter.hh"
#include "Numerics.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include <fstream>
#include <sys/stat.h>

ComptonScatter::ComptonScatter(PhotonField* field, double comMin):
PhotonProcess(field, comMin, "ComptonScatter")
{
    loadCrossSection();
}

#ifdef USEGP
ComptonScatter::ComptonScatter(PhotonField* field, int trainSize,
    double errorMax, double comMin, std::string saveDir):
PhotonProcess(field, comMin, trainSize, errorMax, saveDir, "ComptonScatter")
{
    loadCrossSection();
}

ComptonScatter::ComptonScatter(PhotonField* field, const G4String& gpDir,
    double errorMax, double comMin, std::string saveDir):
PhotonProcess(field, comMin, gpDir, errorMax, saveDir, "ComptonScatter")
{
    loadCrossSection();
}
#endif

ComptonScatter::~ComptonScatter()
{
}

G4bool ComptonScatter::IsApplicable(const G4ParticleDefinition&
    particle) const
{
    return (&particle == G4Positron::Positron()
        ||   &particle == G4Electron::Electron());
}
    
G4VParticleChange* ComptonScatter::PostStepDoIt(const G4Track& aTrack,
    const G4Step& aStep)
{
   /* Get the properties of both the interacting lepton and a sampled 
       photon from the photon field */
    aParticleChange.Initialize(aTrack);
   const G4DynamicParticle *aDynamicLepton = aTrack.GetDynamicParticle(); 
    G4Material* aMaterial = aTrack.GetMaterial();
    int blockID = aMaterial->GetMaterialPropertiesTable()
        ->GetConstProperty("BlockID");

    /* Lepton properties */
    double leptonEnergyIn = aDynamicLepton->GetTotalEnergy();
    G4ThreeVector leptonDirection = aDynamicLepton->GetMomentumDirection();

    /* Rotations to gamma frame */
    G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1).cross(-leptonDirection);
    double rotationAngle = leptonDirection.angle(G4ThreeVector(0, 0, 1));
    m_rotaion  = G4RotationMatrix(rotationAxis, rotationAngle);
    G4RotationMatrix rotateBack = G4RotationMatrix(rotationAxis,
        -rotationAngle);

    /* Photon properties */
    double photonEnergyIn, comEnergy, photonThetaIn, photonPhiIn;
    samplePhotonField(blockID, leptonEnergyIn, photonEnergyIn, comEnergy,
        photonPhiIn);
    photonThetaIn = centreOfMassTheta(comEnergy, leptonEnergyIn,
        photonEnergyIn);

    /* Find particles properties in the COM frame*/
    double photonEnergyOut = electron_mass_c2 * (comEnergy - 1.0)
        / (2.0 * std::sqrt(comEnergy));
    double lepronEnergyOut = electron_mass_c2 * (comEnergy + 1.0)
        / (2.0 * std::sqrt(comEnergy));
    double leptonMomentumOut = std::sqrt(lepronEnergyOut * lepronEnergyOut
        - electron_mass_c2 * electron_mass_c2);
    double photonPhiOut   = 2.0 * pi * G4UniformRand();
    double photonThetaOut = samplePairAngle(comEnergy);

    G4LorentzVector photonVector = G4LorentzVector(
        std::sin(photonThetaOut) * std::cos(photonPhiOut) * photonEnergyOut,
        std::sin(photonThetaOut) * std::sin(photonPhiOut) * photonEnergyOut,
        std::cos(photonThetaOut) * photonEnergyOut, photonEnergyOut);
    G4LorentzVector leptonVector = G4LorentzVector(
        -std::sin(photonThetaOut) * std::cos(photonPhiOut) * leptonMomentumOut,
        -std::sin(photonThetaOut) * std::sin(photonPhiOut) * leptonMomentumOut,
        -std::cos(photonThetaOut) * leptonMomentumOut, lepronEnergyOut);

    double beta = std::sqrt(leptonEnergyIn * leptonEnergyIn -  electron_mass_c2
        * electron_mass_c2 + photonEnergyIn * photonEnergyIn + 2.0
        * std::sqrt(leptonEnergyIn * leptonEnergyIn - electron_mass_c2
            * electron_mass_c2) * photonEnergyIn * std::cos(photonThetaIn))
        / (leptonEnergyIn + photonEnergyIn);

    photonVector.boostZ(beta);
    leptonVector.boostZ(beta);

    /* Apply rotation into frame with lepton along z axis */  
    double leptonTheta = std::atan2(photonEnergyIn * std::sin(photonThetaIn),
        (std::sqrt(leptonEnergyIn * leptonEnergyIn - 1.0) - photonEnergyIn
            * std::cos(photonThetaIn)));
    photonVector.rotateY(leptonTheta);
    leptonVector.rotateY(leptonTheta);

    /* Apply rotation around the gamma axis by -photonPhi */
    photonVector.rotateZ(-photonPhiIn);
    leptonVector.rotateZ(-photonPhiIn);

    /* Apply final rotation into simulation frame */
    photonVector = rotateBack(photonVector);
    leptonVector = rotateBack(leptonVector);

    /* Add new photon and update lepton momentum */
    G4DynamicParticle* photon = new G4DynamicParticle(
            G4Gamma::Gamma(), photonVector);

    // Changes
    aParticleChange.SetNumberOfSecondaries(2);
    G4DynamicParticle* positron = new G4DynamicParticle(
            G4Positron::Positron(), leptonVector);
    aParticleChange.AddSecondary(photon);
    aParticleChange.AddSecondary(positron);
    aParticleChange.ProposeMomentumDirection(G4ThreeVector(0,0,0));
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    /*
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(photon);
    aParticleChange.ProposeMomentumDirection(leptonVector.v().unit());
    aParticleChange.ProposeEnergy(leptonVector[3]);
    */
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

double ComptonScatter::crossSection(double comEnergy) const
{
    return Numerics::interpolate1D(m_comEnergy, m_totalCrossSection,
        comEnergy);
}

double ComptonScatter::centreOfMassEnergy(double dynamicEnergy,
    double staticEnergy, double theta) const
{
    return 1.0 + (2.0 * dynamicEnergy * staticEnergy * (1.0
        - std::sqrt(1.0 - (electron_mass_c2 * electron_mass_c2) 
            / (dynamicEnergy * dynamicEnergy)) * std::cos(theta)))
        / (electron_mass_c2 * electron_mass_c2);
}

double ComptonScatter::centreOfMassStatic(double comEnergy,
    double dynamicEnergy, double theta) const
{ 
    return (comEnergy - 1) * electron_mass_c2 * electron_mass_c2
        / (2.0 * dynamicEnergy) / (1.0 - std::sqrt(1.0 -  electron_mass_c2
            * electron_mass_c2 / (dynamicEnergy * dynamicEnergy)));
}

double ComptonScatter::centreOfMassDynamic(double comEnergy,
    double staticEnergy, double theta) const
{
    return 0;
}

double ComptonScatter::centreOfMassTheta(double comEnergy,
    double dynamicEnergy, double staticEnergy) const
{
    return std::acos((1.0 - (comEnergy - 1) * electron_mass_c2
        * electron_mass_c2 / (2.0 * dynamicEnergy * staticEnergy))
        / std::sqrt(1 - electron_mass_c2 * electron_mass_c2
            / (dynamicEnergy * dynamicEnergy)));
}

void ComptonScatter::loadCrossSection()
{
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