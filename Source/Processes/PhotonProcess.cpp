#include "PhotonProcess.hh"
#include "Numerics.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

PhotonProcess::PhotonProcess(PhotonField* field, double comMin,
    const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field), m_comMin(comMin)
{
}

#ifdef USEGP
PhotonProcess::PhotonProcess(PhotonField* field, double comMin, int trainSize,
    double errorMax, bool save, const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field),
m_comMin(comMin), m_errorMax(errorMax), m_save(save)
{
    int numBlocks = m_field->getNumBlocks();
    int dims = m_field->fieldDimensions();
    m_gp = new GaussianProcess(numBlocks, dims, trainSize);
}

PhotonProcess::PhotonProcess(PhotonField* field,  double comMin,
    const G4String& gpDir, int trainSize, double errorMax, bool save,
    const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field),
m_comMin(comMin), m_errorMax(errorMax), m_save(save)
{
    m_gp = new GaussianProcess(gpDir);
}
#endif

PhotonProcess::~PhotonProcess()
{
#ifdef USEGP  
    delete m_gp;
#endif
}

G4double PhotonProcess::GetMeanFreePath(const G4Track& track, G4double,
         G4ForceCondition*)
{
    /* check if inside a radiation block */
    G4Material* aMaterial = track.GetMaterial();
    if(aMaterial->GetMaterialPropertiesTable()->GetConstProperty("Radiation")
            < 0.0) return 1e99;

    /* get the ID of the block */
    double blockX = aMaterial->GetMaterialPropertiesTable()
            ->GetConstProperty("Xpos");
    double blockY = aMaterial->GetMaterialPropertiesTable()
            ->GetConstProperty("Ypos");
    double blockZ = aMaterial->GetMaterialPropertiesTable()
            ->GetConstProperty("Zpos");
    int position[3] = {(int)blockX, (int)blockY, (int)blockZ};
    int blockID = m_field->getBlockID(position);

    /* Get the interacting particle properties */
    const G4DynamicParticle *aDynamicGamma = track.GetDynamicParticle();
    double dynamicEnergy = aDynamicGamma->GetKineticEnergy();
    G4ThreeVector gammaDirection = aDynamicGamma->GetMomentumDirection();
    double dynamicTheta = std::acos(gammaDirection[2]);
    double dynamicPhi = std::atan2(gammaDirection[1],
        gammaDirection[0] + 1e-99);
    dynamicPhi = dynamicPhi < 0 ? 2 * pi + dynamicPhi : dynamicPhi;

    /* Find the rotation matrices that rotate the gamma ray onto
       the z axis */
    G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1).cross(gammaDirection);
    double rotationAngle = gammaDirection.angle(G4ThreeVector(0, 0, 1));
    m_rotaion = G4RotationMatrix(rotationAxis,
            rotationAngle);

    /* Check if max (head on) COM is too low */
    double comMax = centreOfMassEnergy(dynamicEnergy,
        *m_field->getEnergy().end(), pi);
    if(comMax < m_comMin) return 1e99;

    /* If nonIsotropic, check if field is dense above COM threshold */
    /* Return 1e99 if angle is also too low. */
    if(!m_field->isIsotropic())
    {
        int thetaMinIndex = Numerics::vectorIndex(m_field->getTheta(), thetaMin);
        for (int i = 0; i < angleRes; i++)
        {
            for (int j = thetaMinIndex + 1; j < angleRes; j++)
            {
                double angleIn = {m_field->getTheta()[j],
                    m_field->getPhi()[i]};
                double angleOut[2];
                rotateThetaPhi(angleIn, angleOut, m_rotaion);
                if (Numerics::interpolate2D(m_field->getTheta(),
                    m_field->getPhi(), m_field->getAngleDensity(blockID),
                    angleOut) > m_minDensity)
                {
                   goto isDense;
                }
            }
            if (i == angleRes - 1)
            {
                return 1e99;
            }
        }
    }
    isDense:

#ifdef USEGP
    /* Check if we can use the GP instead */
    double gpOut[2];
    if (m_field->isIsotropic())
    {
        m_gp->run(blockID, {dynamicEnergy}, gpOut);
    } else
    {
        m_gp->run(blockID, {dynamicEnergy, dynamicTheta, dynamicPhi}, gpOut);
    }
    
    if (gpOut[1] < m_errorMax)
    {
        return gpOut[0];
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
                    m_field->getEnergy()[i], m_field->getTheta()[j]);
                if (comEnergy > m_comMin)
                {
                    thetaInt[j] = crossSection(comEnergy) * (1.0 -
                        std::cos(m_field->getTheta()[j]));
                } else 
                {
                    comInt[j] = 0;
                }
            }
            energyInt[i] = m_field->getEnergyDensity()[i]
                * Numerics::simpsons(m_field->getTheta(), thetaInt);
        }
        meanPath = 1.0 / Numerics::simpsons(m_field->getEnergyDensity(),
            energyInt);
    } else // Is anisotropic
    {
        Vector<double> thetaInt(m_field->getAngleRes());
        Vector<double> phiInt(m_field->getAngleRes());
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
                    double angleIn[] = {m_field->getTheta()[j],
                        m_field->getPhi()[k]};
                    double angleOut[2];
                    rotateThetaPhi(angleIn, angleOut, m_rotaion);
                    phiInt[k] = Numerics::interpolate2D(m_field->getTheta(),
                        m_field->getPhi(), m_field->getAngleDensity(blockID),
                        angleOut);
                }
                double comEnergy = centreOfMassEnergy(dynamicEnergy,
                    m_field->getEnergy()[i], m_field->getTheta()[j]);
                if (comEnergy > m_comMin)
                {
                    thetaInt[j] = crossSection(comEnergy) * (1.0 -
                        std::cos(m_field->getTheta()[j]));
                } else 
                {
                    comInt[j] = 0;
                }
            }
            energyInt[i] = m_field->getEnergyDensity()[i]
                * Numerics::simpsons(m_field->getTheta(), thetaInt);       
        }
        meanPath = 1.0 / Numerics::simpsons(m_field->getEnergyDensity(),
            energyInt);
    }
#ifdef USEGP
    if (m_field->isIsotropic())
    {
        m_gp->addData(blockID, {dynamicEnergy}, meanPath);
    } else
    {
        m_gp->addData(blockID, {dynamicEnergy, dynamicTheta, dynamicPhi},
            meanPath);
    }
#endif
    return meanPath;
}

void PhotonProcess::samplePhotonField(int blockID, double dynamicEnergy,
    double& photonEnergy, double& comEnergy, double& photonPhi) const
{
    /* Find mode of distribution to bound sampling. Currently use a slow
       coordinate search here, might be better using a better method */
    /* first for isotrpoic fields */
    double maxDensity(0), currentDensity(0);
    if(m_field->isIsotropic())
    {
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            double s = centreOfMassEnergy(dynamicEnergy,
                m_field->getEnergy()[i], pi);
            if (s > m_comMin)
            {
                currentDensity = 2.0 * crossSection(s)
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
                photonEnergy, pi));
            photonPhi = G4UniformRand() * 2.0 * pi;
            density = 2.0 * crossSection(comEnergy)
                * Numerics::interpolate1D(m_field->getEnergy(),
                    m_field->getEnergyDensity(), photonEnergy);
            randDensity = G4UniformRand() * maxDensity;
        } while (randDensity > density);
    } else
    {
        /* The field is no isotro::pic */
        for (int i = 0; i < m_field->getEnergyRes(); i++) // loop energy
        {
            for (int j = 0; j < m_field->getAngleRes(); j++) // loop theta
            {
                for (int k = 0; k < m_field->getAngleRes(); k++) // loop phi
                {
                    double angleIn[2];
                    double angleOut[2];
                    rotateThetaPhi(angleIn, angleOut, m_rotaion);
                    double s = centreOfMassEnergy(dynamicEnergy,
                        m_field->getEnergy()[i], angleOut[0]);
                    if (s > m_comMin)
                    {
                        currentDensity = crossSection(s)
                            * m_field->getEnergyDensity()[i]
                            * Numerics::interpolate2D(m_field->getTheta(),
                                m_field->getPhi(),
                                m_field->getAngleDensity(blockID), angleOut)
                            * (1.0 - std::cos(angleOut[0]));
                        maxDensity = currentDensity > maxDensity
                            ? currentDensity : maxDensity;
                    }
                }
            }
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
                photonEnergy, pi));
            double photonTheta = centreOfMassTheta(comEnergy, dynamicEnergy,
                photonEnergy);
            photonPhi = m_field->getPhi()[0] + G4UniformRand()
                * (*m_field->getPhi().end() - m_field->getPhi()[0]);
            double angleIn[] = {photonTheta, photonPhi};
            double angleOut[2];
            rotateThetaPhi(angleIn, angleOut, m_rotaion);
            density = crossSection(comEnergy)
                * Numerics::interpolate1D(m_field->getEnergy(),
                    m_field->getEnergyDensity(), photonEnergy)
                * Numerics::interpolate2D(m_field->getTheta(),
                    m_field->getPhi(),
                    m_field->getAngleDensity(blockID), angleOut)
                * (1.0 - std::cos(angleOut[0]));
            randDensity = G4UniformRand() * maxDensity;
        } while (randDensity > density);
    }
}

double PhotonProcess::samplePairAngle(double comEnergy)
{
    double maxDensity = diffCrossSection(comEnergy, pi);
    double randAngle;
    double randDensity;
    double angleDensity;
    do
    {
        randAngle    = std::acos(2.0 * G4UniformRand() - 1.0);
        randDensity  = G4UniformRand() * maxDensity;
        angleDensity = diffCrossSection(comEnergy, randAngle);
    } while (randDensity > angleDensity);
    return randAngle;
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