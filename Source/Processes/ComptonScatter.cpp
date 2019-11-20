#include "ComptonScatter.hh"
#include "Numerics.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"


ComptonScatter::ComptonScatter(PhotonField* field,
    const std::string& dataFile, double comMin):
PhotonProcess(field, comMin, "ComptonScatter")
{
    openDataFile(dataFile);
}

#ifdef USEGP
ComptonScatter::ComptonScatter(PhotonField* field,
    const std::string& dataFile, int trainSize, double errorMax, bool save,
    double comMin):
PhotonProcess(field, comMin, trainSize, errorMax, save, "ComptonScatter")
{
    openDataFile(dataFile);
}

ComptonScatter::ComptonScatter(PhotonField* field,
    const std::string& dataFile, const G4String& gpDir, int trainSize,
    double errorMax, bool save, double comMin):
PhotonProcess(field, comMin, gpDir, trainSize, errorMax,
    save, "ComptonScatter")
{
    openDataFile(dataFile);
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
    double leptonEnergyIn = aDynamicLepton->GetKineticEnergy();
    G4ThreeVector leptonDirection = aDynamicLepton->GetMomentumDirection();

    /* Rotations to gamma frame */
    G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1).cross(leptonDirection);
    double rotationAngle = leptonDirection.angle(G4ThreeVector(0, 0, 1));
    m_rotaion  = G4RotationMatrix(rotationAxis, rotationAngle);
    G4RotationMatrix rotateBack = G4RotationMatrix(rotationAxis,
        -rotationAngle);

    /* Photon properties */
    double photonEnergyIn, comEnergy, photonTheta, photonPhi;
    samplePhotonField(blockID, leptonEnergyIn, photonEnergyIn, comEnergy,
        photonPhi);
    photonTheta = std::acos((((comEnergy - 1.0) * electron_mass_c2
        * electron_mass_c2 ) / (2.0 * leptonEnergyIn * photonEnergyIn) - 1.0)
        / (std::sqrt(1.0 - electron_mass_c2 * electron_mass_c2
            / (leptonEnergyIn * leptonEnergyIn))));

    /* Find particles properties in the COM frame*/
    double photonEnergyOut = (comEnergy - 1.0) / (2.0 * std::sqrt(comEnergy));
    double lepronEnergyOut = (comEnergy + 1.0) / (2.0 * std::sqrt(comEnergy));
    double photonPhiOut   = 2.0 * pi * G4UniformRand();
    double photonThetaOut = samplePairAngle(comEnergy);

    G4LorentzVector photonVector = G4LorentzVector(
        std::sin(photonThetaOut) * std::cos(photonPhiOut),
        std::sin(photonThetaOut) * std::sin(photonPhiOut),
        std::cos(photonThetaOut), photonEnergyOut);
    G4LorentzVector leptonVector = G4LorentzVector(
        -std::sin(photonThetaOut) * std::cos(photonPhiOut),
        -std::sin(photonThetaOut) * std::sin(photonPhiOut),
        -std::cos(photonThetaOut), lepronEnergyOut);

    double beta = std::sqrt((photonEnergyIn * photonEnergyIn
        + leptonEnergyIn * leptonEnergyIn) / (electron_mass_c2
        * electron_mass_c2) - 1.0 - 2.0 * leptonEnergyIn * photonEnergyIn
        * std::sqrt(1.0 - electron_mass_c2 * electron_mass_c2
            / (leptonEnergyIn * leptonEnergyIn))
        * std::cos(photonTheta)) * electron_mass_c2
        / (photonEnergyIn + leptonEnergyIn);
    photonVector.boostZ(-beta);
    leptonVector.boostZ(-beta);

    /* Apply rotation into frame with lepton along z axis */  
    double leptonTheta = std::atan2(photonEnergyIn * std::sin(photonTheta),
        (std::sqrt(leptonEnergyIn * leptonEnergyIn - 1.0) - photonEnergyIn
            * std::cos(photonTheta)));

    photonVector.rotateY(leptonTheta);
    leptonVector.rotateY(leptonTheta);

    /* Apply rotation around the gamma axis by -photonPhi */
    photonVector.rotateZ(-photonPhi);
    leptonVector.rotateZ(-photonPhi);

    /* Apply final rotation into simulation frame */
    photonVector = rotateBack(photonVector);
    leptonVector = rotateBack(leptonVector);

    /* Add new photon and update lepton momentum */
    G4DynamicParticle* photon = new G4DynamicParticle(
            G4Gamma::Gamma(), photonVector);            
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(photon);
    aParticleChange.ProposeMomentumDirection(leptonVector.v().unit());
    aParticleChange.ProposeEnergy(leptonVector[3]);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

double ComptonScatter::crossSection(double comEnergy) const
{
    return Numerics::interpolate1D(m_comEnergy, m_totalCrossSection,
        comEnergy);
}

double ComptonScatter::diffCrossSection(double comEnergy, double theta) const
{
    double queryPoint[2] = {comEnergy, theta};
    return Numerics::interpolate2D(m_comEnergy, m_comTheta,
        m_diffCrossSection, queryPoint);
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

void ComptonScatter::openDataFile(const std::string& fileName)
{
    m_comEnergy.open(fileName, "/ComEnergy");
    m_comTheta.open(fileName, "/Theta");
    m_totalCrossSection.open(fileName, "/TotalCrossSection");
    m_diffCrossSection.open(fileName, "/DiffCrossSection");
}