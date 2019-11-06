#include "PhotonProcess.hh"


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
