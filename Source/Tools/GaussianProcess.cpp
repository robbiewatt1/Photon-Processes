#include "GaussianProcess.hh"
#include <fstream>

GaussianProcess::GaussianProcess(int gpSize, int inputSize, int trainSize):
m_gpSize(gpSize), m_trainSize(trainSize), m_inputSize(inputSize)
{
    m_gausProc = new libgp::GaussianProcess [m_gpSize];
    m_switch = new bool [m_gpSize];
    m_input = new Eigen::VectorXd* [m_gpSize];
    m_output = new double* [m_gpSize];
    
    for (int i = 0; i < m_gpSize; ++i) // loop over GPs
    {
        // allocate memory for training input / output
        m_input[i] = new Eigen::VectorXd[trainSize];
        m_output[i] = new double [trainSize];

        // No GPs have been trained yet
        m_switch[i] = false;

        // Set up the GPs
        libgp::GaussianProcess process
            = libgp::GaussianProcess(m_inputSize, "CovSEiso");
        Eigen::VectorXd params(process.covf().get_param_dim());
        for(int i = 0; i < m_inputSize; i++) params[i] = -1;
        process.covf().set_loghyper(params);
        m_gausProc[i] =  process;
    }
    // Set deult normilisation
    double inputNorm[m_inputSize];
    for (int i = 0; i < m_inputSize; ++i) m_inputNorm[i] = 1.0;
    setNormParams(inputNorm, 1.0);
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
    // Allocate memory
    m_gausProc = new libgp::GaussianProcess [m_gpSize];
    m_switch = new bool [m_gpSize];
    m_input = new Eigen::VectorXd* [m_gpSize];
    m_output = new double* [m_gpSize]

    for (int i = 0; i < m_gpSize; ++i) // loop over GPs
    {
        // allocate memory for training input / output
        m_input[i] = new Eigen::VectorXd[trainSize];
        m_output[i] = new double [trainSize];
        m_gausProc[i] = libgp::GaussianProcess(("./" + gpDir + "/"
            + std::to_string(i) + ".gp").c_str());
    }
    for (int i = 0; i < m_inputSize; ++i) meta >> m_inputNorm[i];
    meta >> m_outputNorm;
    for (int i = 0; i < m_gpSize; ++i) meta >> m_switch[i] << "\t";
}

GaussianProcess::~GaussianProcess()
{
    // dealocate memory
    for (int i = 0; i < m_m_gpSize; ++i)
    {
        delete [] m_input[i];
        delete [] m_output[i];
    }
    delete [] m_gausProc;
    delete [] m_switch;
    delete [] m_input;
    delete [] m_output;
}

void GaussianProcess::save(std::string outputDir)
{
    // First save the meta data
    std::fstream meta(outputDir + "/meta");
    meta << m_gpSize << "\n";
    meta << m_trainSize << "\n";
    meta << inputSize << "\n";
    for (int i = 0; i < m_inputSize; ++i) meta << m_inputNorm[i] << "\t";
    meta << "\n";
    meta << m_outputNorm;
    for (int i = 0; i < m_gpSize; ++i) meta << m_switch[i] << "\t";
    meta.close();

    // Save the GPs
    for (int i = 0; i < m_gpSize; ++i)
    {
        m_gausProc[i].write((outputDir + "/" + std::to_string(i)
            + ".gp").c_str());
    }  
}

void GaussianProcess::setNormParams(double inputNorm[], outputNorm)
{
    m_outputNorm = outputNorm;
    m_inputNorm = double [m_inputSize];
    for (int i = 0; i < m_inputSize; ++i) m_inputNorm[i] = inputNorm[i];
}

void GaussianProcess::run(int gpID, double input[], double output[2])
{
    for (int i = 0; i < m_inputSize; ++i) input[i] /= inputNorm[i];
    output[0] = m_gausProc[gpID].f(input) / m_outputNorm;
    output[1] = m_gausProc[gpId].var(input) / (m_outputNorm * m_outputNorm);
}

void::addData(int gpID, double input[], double output)
{
    Eigen::VectorXd dataPoint(inputSize);
    for (int i = 0; i < inputSize; ++i) dataPoint[i] = input[i];
    m_input[gpID][m_trainCount] = dataPoint;
    m_output[gpID][m_trainCount] = output;
}

void GaussianProcess::train(int gpID)
{
    std::cout << "Optimising GP " << gpID << "..." << std::endl;
    for (int i = 0; i < m_trainSize; i++)
    {
        m_gausProc[gpID].add_pattern(m_input[i], m_output[i]);
    }
    m_optimiser.maximize(&m_gausProc[gpID], 200, 1);
    m_trainCount[gpID] = 0;
    m_switch[gpID] = true;
    std::cout << "Optimisation complete!" << std::endl;
}