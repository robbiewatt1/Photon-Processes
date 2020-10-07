#include "PhotonProcess.hh"
#include "Numerics.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <sys/stat.h>

#include<fstream>

PhotonProcess::PhotonProcess(PhotonField* field, double comMin,
    const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field), m_comMin(comMin),
m_useGP(false), m_multiplier(1)
{
    loadDiffCrossSection();
}

#ifdef USEGP
PhotonProcess::PhotonProcess(PhotonField* field, double comMin, int trainSize,
    double errorMax, std::string saveDir, const G4String& name,
    G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field),
m_comMin(comMin), m_useGP(true), m_multiplier(1), m_errorMax(errorMax),
m_saveDir(saveDir)
{
    int numBlocks = m_field->getNumBlocks();
    int dims = m_field->fieldDimensions();
    m_gp = new GaussianProcess(numBlocks, dims, trainSize);
    loadDiffCrossSection();
}

PhotonProcess::PhotonProcess(PhotonField* field,  double comMin,
    const G4String& gpDir, double errorMax, std::string saveDir,
    const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field),
m_comMin(comMin), m_useGP(true), m_multiplier(1), m_errorMax(errorMax),
m_saveDir(saveDir)
{
    m_gp = new GaussianProcess(gpDir);
    loadDiffCrossSection();
}
#endif

PhotonProcess::~PhotonProcess()
{
#ifdef USEGP
    if(m_useGP)
    {
        if (m_saveDir != "")
        {
            m_gp->save(m_saveDir);
        }
        delete m_gp;
    }
#endif
}

#ifdef USEGP
void PhotonProcess::setParamsGP(const Vector<double>& inputNorm,
    double outputNorm)
{
    m_gp->setNormParams(inputNorm, outputNorm);
}
#endif

G4double PhotonProcess::GetMeanFreePath(const G4Track& track, G4double,
         G4ForceCondition*)
{
    /* check if inside a radiation block */
    G4Material* aMaterial = track.GetMaterial();
    if(aMaterial->GetMaterialPropertiesTable()->GetConstProperty("Radiation")
            < 0.0) return 1e99;
    int blockID = aMaterial->GetMaterialPropertiesTable()
        ->GetConstProperty("BlockID");

    /* Get the interacting particle properties */
    const G4DynamicParticle *dynamicParticle = track.GetDynamicParticle();
    double dynamicEnergy = dynamicParticle->GetKineticEnergy();
    G4ThreeVector particleDirection = dynamicParticle->GetMomentumDirection();
    double dynamicTheta = std::acos(particleDirection[2]);
    double dynamicPhi = std::atan2(particleDirection[1],
        particleDirection[0] + 1e-99);
    dynamicPhi = dynamicPhi < 0 ? 2 * pi + dynamicPhi : dynamicPhi;

    /* Find the rotation matrices that rotate the gamma ray onto
       the z axis */
    G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1)
        .cross(-particleDirection);
    double rotationAngle = particleDirection.angle(G4ThreeVector(0, 0, 1));
    m_rotaion = G4RotationMatrix(rotationAxis,
            rotationAngle);
    /* Check if max (head on) COM is too low */
    double comMax = centreOfMassEnergy(dynamicEnergy,
        *m_field->getEnergy().end(), *m_field->getTheta(blockID).end());
    if(comMax < m_comMin) return 1e99;

#ifdef USEGP
    if (m_useGP)
    {
        /* Check if we can use the GP instead */
        double gpOut[2];
        if (m_field->isIsotropic())
        {
            m_gp->run(blockID, {dynamicEnergy}, gpOut);
        } else
        {
            m_gp->run(blockID, {dynamicEnergy, dynamicTheta, dynamicPhi},
                gpOut);
        }
        if (std::sqrt(gpOut[1]) / gpOut[0] < m_errorMax)
        {
            return gpOut[0] / m_multiplier;
        }
    }
#endif

    /* Long integration */
    double meanPath;
    if (m_field->isIsotropic())
    {
        Vector<double> thetaInt(m_field->getAngleRes());
        Vector<double> energyInt(m_field->getEnergyRes());
        // Integral over photon energy epsilon
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            // Integrate over theta
            for (int j = 0; j < m_field->getAngleRes(); j++)
            {
                double comEnergy = centreOfMassEnergy(dynamicEnergy,
                    m_field->getEnergy()[i], m_field->getTheta(blockID)[j]);
                if (comEnergy > m_comMin)
                {
                    thetaInt[j] = crossSection(comEnergy)
                        * (1.0 - std::cos(m_field->getTheta(blockID)[j]))
                        * std::sin(m_field->getTheta(blockID)[j]) / 2.0;
                } else
                {
                    thetaInt[j] = 0;
                }
            }
            energyInt[i] = m_field->getEnergyDensity()[i]
                * Numerics::simpsons(m_field->getTheta(blockID), thetaInt);
        }
        meanPath = 1.0 / (classic_electr_radius * classic_electr_radius
            * Numerics::simpsons(m_field->getEnergy(), energyInt));
    } else // Is anisotropic
    {
        Vector<double> phiInt(m_field->getAngleRes());
        Vector<double> thetaInt(m_field->getAngleRes());
        Vector<double> energyInt(m_field->getEnergyRes());
        // Integral over photon energy epsilon
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            // Integrate over theta
            for (int j = 0; j < m_field->getAngleRes(); j++)
            {
                // Integrate over phi
                for (int k = 0; k < m_field->getAngleRes(); k++)
                {
                    double angleIn[] = {m_field->getTheta(blockID)[j],
                        m_field->getPhi(blockID)[k]};
                    double angleOut[2];
                    rotateThetaPhi(angleIn, angleOut, m_rotaion);
                    double comEnergy = centreOfMassEnergy(dynamicEnergy,
                        m_field->getEnergy()[i], angleOut[0]);
                    if (comEnergy > m_comMin)
                    {
                        phiInt[k] = m_field->getAngleDensity(blockID)[j][k]
                            * crossSection(comEnergy)
                            * (1.0 - std::cos(angleOut[0]));
                    } else
                    {
                        phiInt[k] = 0;
                    }
                }
                thetaInt[j] = Numerics::simpsons(m_field->getPhi(blockID), phiInt);
            }
            energyInt[i] = m_field->getEnergyDensity()[i]
                * Numerics::simpsons(m_field->getTheta(blockID), thetaInt);
        }
        meanPath = 1.0 / (classic_electr_radius * classic_electr_radius
            * Numerics::simpsons(m_field->getEnergy(), energyInt));
    }
#ifdef USEGP
    if (m_useGP)
    {
        if (m_field->isIsotropic())
        {
            m_gp->addData(blockID, {dynamicEnergy}, meanPath);
        } else
        {
            m_gp->addData(blockID, {dynamicEnergy, dynamicTheta, dynamicPhi},
                meanPath);
        }
    }
#endif
    //std::cout << meanPath / m_multiplier << std::endl;
    return meanPath / m_multiplier;
}

void PhotonProcess::samplePhotonField(int blockID, double dynamicEnergy,
    double& photonEnergy, double& comEnergy, double& photonPhi) const
{
    /* Find mode of distribution to bound sampling. Currently use a slow
       coordinate search here, might be better using a better method */
    /* first for isotrpoic fields */
    if(m_field->isIsotropic())
    {
        double maxDensity(0), currentDensity(0);
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            double s = centreOfMassEnergy(dynamicEnergy,
                m_field->getEnergy()[i], pi);
            if (s > m_comMin)
            {
                currentDensity = crossSection(s)
                    * m_field->getEnergyDensity()[i]; 
            }
            maxDensity = currentDensity > maxDensity ? currentDensity
                        : maxDensity;
        }
        double photonEnergyMin = centreOfMassStatic(m_comMin,
            dynamicEnergy, pi);
        double comEnergyMax = centreOfMassEnergy(dynamicEnergy,
            *m_field->getEnergy().end(), pi);
        double density, randDensity;
        do
        {   // While random density is too large
            do
            {   // While s, e combination not possible
                photonEnergy = photonEnergyMin + G4UniformRand()
                        * (*m_field->getEnergy().end() - photonEnergyMin);
                comEnergy = m_comMin + G4UniformRand()
                    * (comEnergyMax - m_comMin);
            } while (comEnergy > centreOfMassEnergy(dynamicEnergy,
                photonEnergy, pi) || comEnergy <
                centreOfMassEnergy(dynamicEnergy, photonEnergy, 0));
            photonPhi = G4UniformRand() * 2.0 * pi;
            density = crossSection(comEnergy)
                * Numerics::interpolate1D(m_field->getEnergy(),
                    m_field->getEnergyDensity(), photonEnergy);
            randDensity = G4UniformRand() * maxDensity;
        } while (randDensity > density);
    } else
    {
        /* The field is no isotropic */
        double maxDensity(0), currentDensity(0), maxTheta(0), minTheta(pi),
            maxPhi(0), minPhi(2 * pi);
        for (int i = 0; i < m_field->getEnergyRes(); i++) // loop energy
        {
            for (int j = 0; j < m_field->getAngleRes(); j++) // loop theta
            {
                for (int k = 0; k < m_field->getAngleRes(); k++) // loop phi
                {
                    double angleIn[] = {m_field->getTheta(blockID)[j],
                        m_field->getPhi(blockID)[k]};
                    double angleOut[2];
                    rotateThetaPhi(angleIn, angleOut, m_rotaion);
                    double comEnergy = centreOfMassEnergy(dynamicEnergy,
                        m_field->getEnergy()[i], angleOut[0]);
                    currentDensity = m_field->getAngleDensity(blockID)[j][k]
                        * crossSection(comEnergy)
                        * (1.0 - std::cos(angleOut[0]))
                        * m_field->getEnergyDensity()[i];
                    maxTheta = angleOut[0] > maxTheta
                        ? angleOut[0] : maxTheta;
                    minTheta = angleOut[0] < minTheta
                        ? angleOut[0] : minTheta;
                    maxPhi = angleOut[1] > maxPhi
                        ? angleOut[1] : maxPhi;
                    minPhi = angleOut[1] < minPhi
                        ? angleOut[1] : minPhi;
                    maxDensity = currentDensity > maxDensity
                        ? currentDensity : maxDensity;
                }
            }
        }
        // Find the limits for sampling from
        double photonEnergyMin = centreOfMassStatic(m_comMin,
            dynamicEnergy, maxTheta);
        photonEnergyMin = m_field->getEnergy()[0] > photonEnergyMin
            ? m_field->getEnergy()[0] : photonEnergyMin;
        double comEnergyMin = centreOfMassEnergy(dynamicEnergy,
            m_field->getEnergy()[0], minTheta);
        comEnergyMin = m_comMin > comEnergyMin
            ? m_comMin : comEnergyMin;
        double comEnergyMax = centreOfMassEnergy(dynamicEnergy,
            *m_field->getEnergy().end(), maxTheta);
        double density, randDensity;

        do
        {   // While random density is too large
            do
            {   // While s, e combination not possible
                photonEnergy = photonEnergyMin + G4UniformRand()
                        * (*m_field->getEnergy().end() - photonEnergyMin);

                // Check if point/laser 
                if (m_field->getAngleRes() == 1)
                {
                    comEnergy = centreOfMassEnergy(dynamicEnergy, photonEnergy,
                        maxTheta);
                } else
                {
                    comEnergy = comEnergyMin + G4UniformRand()
                        * (comEnergyMax - comEnergyMin);
                }
            } while (comEnergy > centreOfMassEnergy(dynamicEnergy,
                photonEnergy, maxTheta) || comEnergy
                    < centreOfMassEnergy(dynamicEnergy, photonEnergy,
                        minTheta));
            double photonTheta = centreOfMassTheta(comEnergy, dynamicEnergy,
                photonEnergy);
            photonPhi = minPhi + G4UniformRand() * (maxPhi - minPhi);
            double angleIn[] = {photonTheta, photonPhi};
            double angleOut[2];
            rotateThetaPhi(angleIn, angleOut, m_rotaion);
            density = crossSection(comEnergy)
                * Numerics::interpolate1D(m_field->getEnergy(),
                    m_field->getEnergyDensity(), photonEnergy)
                * Numerics::interpolate2D(m_field->getTheta(blockID),
                    m_field->getPhi(blockID),
                    m_field->getAngleDensity(blockID), angleOut)
                * (1.0 - std::cos(angleOut[0]));
            randDensity = G4UniformRand() * maxDensity;
        } while (randDensity > density);
    }
}

double PhotonProcess::samplePairAngle(double comEnergy)
{
    if (comEnergy > *m_invCdfTable.xAxis.end())
    {
        std::cerr << "Warning: s is outside of differential cross-section "
                      "table limits for " << theProcessName
                  << ". \nThis may cause a crash!" << std::endl;
    }
    double queryPoint[] = {comEnergy, G4UniformRand()};
    return Numerics::interpolate2D(m_invCdfTable.xAxis, m_invCdfTable.yAxis,
        m_invCdfTable.data, queryPoint);
}

void PhotonProcess::rotateThetaPhi(double angleIn[2], double angleOut[2],
        const G4RotationMatrix& rotaion) const
{
    G4ThreeVector vector = G4ThreeVector(std::sin(angleIn[0])
        * std::cos(angleIn[1]), std::sin(angleIn[0]) * std::sin(angleIn[1]),
        std::cos(angleIn[0]));
    G4ThreeVector vectorPrime = rotaion(vector);
    angleOut[0] = std::acos(vectorPrime[2]);
    angleOut[1] = std::atan2(vectorPrime[1], (vectorPrime[0] + 1e-99));
}

void PhotonProcess::loadDiffCrossSection()
{
    std::string fileName = theProcessName + "_diff.h5";
    // Check if new file exisits in current dir
    struct stat buffer;
    if (stat (fileName.c_str(), &buffer) == 0)
    {
        // Load file from run dir
        m_invCdfTable.xAxis.open(fileName, "/s");
        m_invCdfTable.yAxis.open(fileName, "/p");
        m_invCdfTable.data.open(fileName, "/angle");
    } else
    {
        // Load file from Install
        try
        {
            std::string dataDir(getenv("Photon_Process_Data_Dir"));
            std::string filePath = dataDir + "/" + fileName;
            m_invCdfTable.xAxis.open(filePath, "/s");
            m_invCdfTable.yAxis.open(filePath, "/p");
            m_invCdfTable.data.open(filePath, "/angle");
        } catch (const std::exception& e)
        {
            std::cout << e.what() << std::endl;
            std::cout << "You probably haven't run PhotonProcess.sh"
                << std::endl;
            std::abort();
        }
    }
}
