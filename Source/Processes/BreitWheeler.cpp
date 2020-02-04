#include "BreitWheeler.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

BreitWheeler::BreitWheeler(PhotonField* field, double comMin):
PhotonProcess(field, comMin, "BreitWheeler")
{
}

#ifdef USEGP
BreitWheeler::BreitWheeler(PhotonField* field, int trainSize,
    double errorMax, double comMin, std::string saveDir):
PhotonProcess(field, comMin, trainSize, errorMax, saveDir, "BreitWheeler")
{
}

BreitWheeler::BreitWheeler(PhotonField* field, const G4String& gpDir,
    double errorMax, double comMin, std::string saveDir):
PhotonProcess(field, comMin, gpDir, errorMax, saveDir,
    "BreitWheeler")
{
}
#endif

BreitWheeler::~BreitWheeler()
{
}

G4bool BreitWheeler::IsApplicable(const G4ParticleDefinition&
        particle) const
{
    return (&particle == G4Gamma::Gamma());
}

G4VParticleChange* BreitWheeler::PostStepDoIt(const G4Track& aTrack,
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
    double pairEnergy = electron_mass_c2 * std::sqrt(comEnergy);
    double pairMomentum = std::sqrt(pairEnergy * pairEnergy
        - electron_mass_c2 * electron_mass_c2);
    double pairPhi   = 2.0 * pi * G4UniformRand();
    double pairTheta = samplePairAngle(comEnergy);

    G4LorentzVector electronVector = G4LorentzVector(
        std::sin(pairTheta) * std::cos(pairPhi) * pairMomentum,
        std::sin(pairTheta) * std::sin(pairPhi) * pairMomentum,
        std::cos(pairTheta) * pairMomentum, pairEnergy);
    G4LorentzVector positronVector = G4LorentzVector(
        -std::sin(pairTheta) * std::cos(pairPhi) * pairMomentum,
        -std::sin(pairTheta) * std::sin(pairPhi) * pairMomentum,
        -std::cos(pairTheta) * pairMomentum, pairEnergy);

    double beta = std::sqrt((gammaEnergy * gammaEnergy + photonEnergy
        * photonEnergy + 2.0 * gammaEnergy * photonEnergy
        * std::cos(photonTheta)) / (gammaEnergy * gammaEnergy + photonEnergy
            * photonEnergy + 2.0 * gammaEnergy * photonEnergy));
    electronVector.boostZ(beta);
    positronVector.boostZ(beta);

    /* Apply rotation into frame with gamma along z axis */  
    double thetaGamma = std::atan2(photonEnergy * std::sin(photonTheta),
                 gammaEnergy + photonEnergy * std::cos(photonTheta));
    electronVector.rotateY(thetaGamma);
    positronVector.rotateY(thetaGamma);

    /* Apply rotation around the gamma axis by -photonPhi */
    electronVector.rotateZ(-photonPhi);
    positronVector.rotateZ(-photonPhi);

    /* Apply final rotation into simulation frame */
    electronVector = rotateBack(electronVector);
    positronVector = rotateBack(positronVector);

    /* Add new electron and positron and kill the gamma ray */
    G4DynamicParticle* electron = new G4DynamicParticle(
            G4Electron::Electron(), electronVector);            
    G4DynamicParticle* positron = new G4DynamicParticle(
            G4Positron::Positron(), positronVector);
    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(electron);
    aParticleChange.AddSecondary(positron);
    aParticleChange.ProposeMomentumDirection(G4ThreeVector(0,0,0));
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

double BreitWheeler::crossSection(double comEnergy) const
{
    if (comEnergy < 1.0)
    {
        return 0;
    } else
    {
        double beta = std::sqrt(1.0 - 1.0 / comEnergy);
        return (1.0 - beta * beta) * ((3.0 - beta * beta * beta * beta)
                * std::log((1.0 + beta) / (1.0 - beta)) - 2.0 * beta 
                * (2.0 - beta * beta));
    }
}

double BreitWheeler::diffCrossSection(double comEnergy, double theta) const
{
    double beta = std::sqrt(1.0 - 1.0 / comEnergy);
    double sinT = std::sin(theta);
    double cosT = std::cos(theta);
    return (beta / comEnergy) * (1.0 + 2.0 * beta * beta * sinT * sinT
            - beta * beta * beta * beta - beta * beta * beta * beta
            * sinT * sinT * sinT * sinT) / ((1.0 - beta * beta * cosT * cosT)
            * (1.0 - beta * beta * cosT * cosT));
}

double BreitWheeler::centreOfMassEnergy(double dynamicEnergy,
    double staticEnergy, double theta) const
{
    return dynamicEnergy * staticEnergy * (1.0 - std::cos(theta))
        / (2.0 * electron_mass_c2 * electron_mass_c2);
}

double BreitWheeler::centreOfMassStatic(double comEnergy,
    double dynamicEnergy, double theta) const
{
    return 2.0 * electron_mass_c2 * electron_mass_c2 * comEnergy
        / (dynamicEnergy * (1.0 - std::cos(theta)));
}

double BreitWheeler::centreOfMassDynamic(double comEnergy,
        double staticEnergy, double theta) const
{
    return 2.0 * electron_mass_c2 * electron_mass_c2 * comEnergy
        / (staticEnergy * (1.0 - std::cos(theta)));
}

double BreitWheeler::centreOfMassTheta(double comEnergy,
    double dynamicEnergy, double staticEnergy) const
{
    return std::acos(1.0 - 2.0 * electron_mass_c2 * electron_mass_c2
        * comEnergy / (dynamicEnergy * staticEnergy));
}