#include "GaussianProcess.hh"
#include <fstream>

GaussianProcess::GaussianProcess(int gpSize, int inputSize, int trainSize):
m_gpSize(gpSize), m_trainSize(trainSize), m_inputSize(inputSize)
{
    m_gausProc = Vector<libgp::GaussianProcess*>(m_gpSize);
    m_switch = Vector<bool>(m_gpSize);
    m_input = Matrix<Eigen::VectorXd>(m_gpSize, m_trainSize);
    m_output = Matrix<double>(m_gpSize, m_trainSize);
    
    for (int i = 0; i < m_gpSize; ++i) // loop over GPs
    {
        // No GPs have been trained yet
        m_switch[i] = false;
        // Set up the GPs
        libgp::GaussianProcess* process
            = new libgp::GaussianProcess(m_inputSize, "CovSEiso");
        Eigen::VectorXd params(process->covf().get_param_dim());
        for(int i = 0; i < m_inputSize; i++) params[i] = -1;
        process->covf().set_loghyper(params);
        m_gausProc[i] = process;
    }
    // Set deult normilisation
    Vector<double> inputNorm = Vector<double>(m_inputSize);
    for (int i = 0; i < m_inputSize; ++i) inputNorm[i] = 1.0;
    setNormParams(inputNorm, 1.0);
    m_optimiser.init();
}

GaussianProcess::GaussianProcess(std::string gpDir)
{

    // Open meta file
    std::ifstream meta(gpDir + "/meta");
    if (!meta)
    {
        // Cant find meta
        std::cerr << "Error: Cannot open meta file \"" << gpDir + "/meta"
            << "\".";
        exit(1);
    }
    // load sizes
    meta >> m_gpSize;
    meta >> m_trainSize;
    meta >> m_inputSize;
    m_gausProc = Vector<libgp::GaussianProcess*>(m_gpSize);
    m_switch = Vector<bool>(m_gpSize);
    m_input = Matrix<Eigen::VectorXd>(m_gpSize, m_trainSize);
    m_output = Matrix<double>(m_gpSize, m_trainSize);

    for (int i = 0; i < m_gpSize; ++i) // loop over GPs
    {
        // allocate memory for training input / output
        m_gausProc[i] = new libgp::GaussianProcess(("./" + gpDir + "/"
            + std::to_string(i) + ".gp").c_str());
    }
    for (int i = 0; i < m_inputSize; ++i) meta >> m_inputNorm[i];
    meta >> m_outputNorm;
    for (int i = 0; i < m_gpSize; ++i) meta >> m_switch[i];
}

GaussianProcess::~GaussianProcess()
{
    for (int i = 0; i < m_gpSize; ++i)
    {
        delete m_gausProc[i];
    }
}

void GaussianProcess::save(std::string outputDir)
{
    // First save the meta data
    std::fstream meta(outputDir + "/meta");
    meta << m_gpSize << "\n";
    meta << m_trainSize << "\n";
    meta << m_inputSize << "\n";
    for (int i = 0; i < m_inputSize; ++i) meta << m_inputNorm[i] << "\t";
    meta << "\n";
    meta << m_outputNorm;
    for (int i = 0; i < m_gpSize; ++i) meta << m_switch[i] << "\t";
    meta.close();

    // Save the GPs
    for (int i = 0; i < m_gpSize; ++i)
    {
        m_gausProc[i]->write((outputDir + "/" + std::to_string(i)
            + ".gp").c_str());
    }  
}

void GaussianProcess::setNormParams(const Vector<double>& inputNorm,
    double outputNorm)
{
    m_outputNorm = outputNorm;
    m_inputNorm = Vector<double>(m_inputSize);
    for (int i = 0; i < m_inputSize; ++i) m_inputNorm[i] = inputNorm[i];
}

void GaussianProcess::run(int gpID, const Vector<double>& input,
    double output[2])
{
    if (!m_switch[gpID])
    {   // Not trained yet
        output[0] = 0;
        output[1] = 1e99;
    } else
    {
        for (int i = 0; i < m_inputSize; ++i) input[i] /= m_inputNorm[i];
        output[0] = m_gausProc[gpID]->f(input.begin()) / m_outputNorm;
        output[1] = m_gausProc[gpID]->var(input.begin())
            / (m_outputNorm * m_outputNorm);
    }
}

void GaussianProcess::addData(int gpID, const Vector<double>& input,
    double output)
{
    Eigen::VectorXd dataPoint(m_inputSize);
    for (int i = 0; i < m_inputSize; ++i) dataPoint[i] = input[i];
    m_input[gpID][m_trainCount[gpID]] = dataPoint;
    m_output[gpID][m_trainCount[gpID]] = output;
}

void GaussianProcess::train(int gpID)
{
    std::cout << "Optimising GP " << gpID << "..." << std::endl;
    for (int i = 0; i < m_trainSize; i++)
    {
        m_gausProc[gpID]->add_pattern(m_input[gpID][i], m_output[gpID][i]);
    }
    m_optimiser.maximize(m_gausProc[gpID], 200, 1);
    m_trainCount[gpID] = 0;
    m_switch[gpID] = true;
    std::cout << "Optimisation complete!" << std::endl;
}