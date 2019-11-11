#include "PhotonProcess.hh"
#include "Numerics.hh"

PhotonProcess::PhotonProcess(PhotonField* field, double comMin,
    const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field), m_comMin(comMin)
{
}

#ifdef USEGP
PhotonProcess::PhotonProcess(PhotonField* field, double comMin, int trainSize,
    double errorMax, bool save, const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field), m_errorMax(errorMax),
m_comMin(comMin), m_save(save)
{
    int numBlocks = m_field->getNumBlocks();
    int dims = m_field->fieldDimensions();
    m_gp = new GaussianProcess(numBlocks, dims, trainSize);
}

PhotonProcess::PhotonProcess(PhotonField* field,  double comMin,
    const G4String& gpDir, int trainSize, double errorMax, bool save,
    const G4String& name, G4ProcessType type):
G4VDiscreteProcess(name, type), m_field(field), m_errorMax(errorMax),
m_comMin(comMin), m_save(save)
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

G4double PhotonProcess::PhotonProcess(const G4Track& track, G4double,
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
    double gammaEnergy = aDynamicGamma->GetKineticEnergy();
    G4ThreeVector gammaDirection = aDynamicGamma->GetMomentumDirection();
    double gammaTheta = std::acos(gammaDirection[2]);
    double gammaPhi = std::atan2(gammaDirection[1], gammaDirection[0] + 1e-99);
    gammaPhi = gammaPhi < 0 ? 2 * pi + gammaPhi : gammaPhi;

    /* Find the rotation matrices that rotate the gamma ray onto
       the z axis */
    G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1).cross(gammaDirection);
    double rotationAngle = gammaDirection.angle(G4ThreeVector(0, 0, 1));
    rotateForward = G4RotationMatrix(rotationAxis,
            rotationAngle);

    /* Check if max (head on) COM is too low */
    double comMax = *m_field->getEnergy().end() *  gammaEnergy
        / (electron_mass_c2 * electron_mass_c2);
    if(comMax < m_comMin) return 1e99;

    /* If nonIsotropic, check if field is dense above COM threshold */
    /* Return 1e99 if angle is also too low. */
    if(!m_field->isIsotropic())
    {
        int thetaMinIndex = arrayIndex(photonTheta, thetaMin, angleRes);
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
        m_gp->run(blockID, {gammaEnergy}, gpOut);
    } else
    {
        m_gp->run(blockID, {gammaEnergy, gammaTheta, gammaPhi}, gpOut);
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
        Vector<double> comInt(m_field->getAngleRes());
        Vector<double> comAxis(m_field->getAngleRes());
        Vector<double> energyInt(m_field->getEnergyRes());
        // Integral over photon energy epsilon
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            // Integrate over s
            for (int j = 0; j < m_field->getAngleRes(); j++)
            {
                comAxis[j] = m_field->getEnergy()[i] * gammaEnergy * (1.0 -
                        std::cos(m_field->getEnergy()[i]))
                    / (2.0 * electron_mass_c2 * electron_mass_c2);
                if (comAxis[j] > m_comMin)
                {
                    comInt[j] = crossSection(comAxis[j]) * comAxis[j];
                } else 
                {
                    comInt[j] = 0;
                }
            }

            if (m_field->getEnergy()[i] > electron_mass_c2 * electron_mass_c2
                        / gammaEnergy)
            {
                energyInt[i] = m_field->getEnergyDensity()[i]
                    * Numerics::simpsons(comAxis, comInt) 
                    / (m_field->getEnergy()[i] * m_field->getEnergy()[i]);
            } else
            {
                energyInt[i] = 0;
            }
        }
        meanPath = gammaEnergy * gammaEnergy
            / (Numerics::simpsons(m_field->getEnergyDensity(), energyInt,
                energyRes) * electron_mass_c2 * electron_mass_c2
                * electron_mass_c2 * electron_mass_c2 * classic_electr_radius
                * classic_electr_radius * pi);
    } else // Is anisotropic
    {
        Vector<double> comInt(m_field->getAngleRes());
        Vector<double> comAxis(m_field->getAngleRes());
        Vector<double> energyInt(m_field->getEnergyRes());
        Vector<double> energyInt(m_field->getEnergyRes());

        // Integral over photon energy epsilon
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            // Integrate over s
            for (int j = 0; j < m_field->getAngleRes(); j++)
            {
                for (int k = 0; k < m_field->getAngleRes(); k++)
                {

                    double angleIn[] = {m_field->getTheta()[j],
                        m_field->getPhi()[k]};
                    double angleOut[2];
                    thetaInt[k] = Numerics::interpolate2D(m_field->getTheta(),
                        m_field->getPhi(), m_field->getAngleDensity(blockID),
                        angleOut);
                }

                comAxis[j] = m_field->getEnergy()[i] * gammaEnergy * (1.0 -
                            std::cos(m_field->getEnergy()[i]))
                        / (2.0 * electron_mass_c2 * electron_mass_c2);
                if (comAxis[j] > m_comMin)
                {
                    comInt[j] = crossSection(comAxis[j]) * comAxis[j]
                            * simpsons(m_field->getPhi(), phiInt);
                } else 
                {
                    comInt[j] = 0;
                }
            }

            if (m_field->getEnergy()[i] > electron_mass_c2 * electron_mass_c2
                        / gammaEnergy)
            {
                energyInt[i] = m_field->getEnergyDensity()[i]
                    * Numerics::simpsons(comAxis, comInt) 
                    / (m_field->getEnergy()[i] * m_field->getEnergy()[i]);
            } else
            {
                energyInt[i] = 0;
            }
        }
        meanPath = gammaEnergy * gammaEnergy
            / (Numerics::simpsons(m_field->getEnergyDensity(), energyInt,
                energyRes) * electron_mass_c2 * electron_mass_c2
                * electron_mass_c2 * electron_mass_c2 * classic_electr_radius
                * classic_electr_radius * pi);
    }
#ifdef USEGP
    if (m_field->isIsotropic())
    {
        m_gp->addData(blockID, {gammaEnergy}, meanPath);
    } else
    {
        m_gp->run(blockID, {gammaEnergy, gammaTheta, gammaPhi}, meanPath);
    }
#endif
    return meanPath;
}




void PhotonProcess::SamplePhotonField(int blockID, double gammaEnergy,
    double& photonEnergy, double& comEnergy, double& photonPhi)
{
    /* Find mode of distribution to bound sampling. Currently use a slow
       coordinate search here, might be better using a better method */
    /* first for isotrpoic fields */
    if(m_field->isIsotropic())
    {
        for (int i = 0; i < m_field->getEnergyRes(); i++)
        {
            double s = gammaEnergy * m_field->getEnergy()[i] /
                (electron_mass_c2 * electron_mass_c2);
            currentDensity = s * crossSection(s)
                * m_field->getEnergyDensity()[i]
                / (m_field->getEnergy()[i] * m_field->getEnergy()[i]);
            maxDensity = currentDensity > maxDensity ? currentDensity
                        : maxDensity;
        }
        double photonEnergyMin = electron_mass_c2 * electron_mass_c2 /
            gammaEnergy;
        double comEnergyMax = *m_field->getEnergy().end() * gammaEnergy
            / (electron_mass_c2 * electron_mass_c2);
        double density, randDensity;
        do
        {   // While random density is too large
            do
            {   // While s, e combination not possible
                photonEnergy = photonEnergyMin + G4UniformRand()
                        * (*m_field->getEnergy().end() - photonEnergyMin);
                comEnergy = m_comMin + G4UniformRand()
                    * (comEnergyMax - m_comMin);
            } while (comEnergy > photonEnergy * gammaEnergy
                / (electron_mass_c2 * electron_mass_c2));
            photonPhi = G4UniformRand() * 2.0 * pi;
            density = comEnergy * crossSection(comEnergy)
                * Numerics::interpolate1D(m_field->getEnergy(),
                    m_field->getEnergyDensity(), photonEnergy)
                / (photonEnergy * photonEnergy * maxDensity)
            randDensity = G4UniformRand();
        } while (randDensity > density);
    } else
    {
        /* The field is no isotropic */
        for (int i = 0; i < m_field->getEnergyRes(); i++) // loop energy
        {
            for (int j = 0; j < m_field->getAngleRes(); j++) // loop theta
            {
                for (int k = 0; k < m_field->getAngleRes(); k++) // loop phi
                {
                    double angleIn[2];
                    double angleOut[2];
                    rotateThetaPhi(angleIn, angleOut, m_rotaion);
                    double s = gammaEnergy * m_field->getEnergy()[i]
                        * (1.0 - std::cos(angleOut[0]))
                            / (2.0 * electron_mass_c2 * electron_mass_c2);     
                    if (s > m_comMin)
                    {
                        currentDensity = s * crossSection(s)
                            * Numerics::interpolate1D(m_field->getEnergy(),
                                m_field->getEnergyDensity(), photonEnergy)
                            * Numerics::interpolate2D(m_field->getTheta(),
                                m_field->getPhi(),
                                m_field->getAngleDensity(blockID), angleOut)
                            / (m_field->getEnergy()[i]
                                * m_field->getEnergy()[i]);
                        maxDensity = currentDensity > maxDensity
                            ? currentDensity : maxDensity;
                    }
                }
            }
        }
        double photonEnergyMin = electron_mass_c2 * electron_mass_c2 /
            gammaEnergy;
        double comEnergyMax = *m_field->getEnergy().end() * gammaEnergy
            / (electron_mass_c2 * electron_mass_c2);
        double density, randDensity;
        do
        {   // While random density is too large
            do
            {   // While s, e combination not possible
                photonEnergy = photonEnergyMin + G4UniformRand()
                        * (*m_field->getEnergy().end() - photonEnergyMin);
                comEnergy = m_comMin + G4UniformRand()
                    * (comEnergyMax - m_comMin);
            } while (comEnergy > photonEnergy * gammaEnergy
                / (electron_mass_c2 * electron_mass_c2));
            photonTheta = std::acos(1.0 - 2.0 * electron_mass_c2 * electron_mass_c2
                * comEnergy / (gammaEnergy * photonEnergy));
            photonPhi = m_field->getPhi()[0] + G4UniformRand()
                * (*m_field->getPhi().end - m_field->getPhi()[0]);
            density = comEnergy * crossSection(comEnergy)
                * Numerics::interpolate1D(m_field->getEnergy(),
                    m_field->getEnergyDensity(), photonEnergy)
                * Numerics::interpolate2D(m_field->getTheta(),
                    m_field->getPhi(),
                    m_field->getAngleDensity(blockID),
                    {photonTheta, photonPhi})
                / (photonEnergy * photonEnergy * maxDensity);
            randDensity = G4UniformRand();
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
        const G4RotationMatrix& rotaion);
{
    G4ThreeVector vector = G4ThreeVector(std::sin(angleIn[0])
        * std::cos(angleIn[1]), std::sin(angleIn[0]) * std::sin(angleIn[1]),
        std::cos(angleIn[0]));
    G4ThreeVector vectorPrime = rotaion(vector);
    angleOut[0] = std::acos(vectorPrime[2]);
    angleOut[1] = std::atan2(vectorPrime[1], (vectorPrime[0] + 1e-99));
}