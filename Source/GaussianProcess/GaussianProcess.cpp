#include "GaussianProcess.hh"
#include <fstream>
#include <cmath>

GaussianProcess::GaussianProcess(int gpSize, int inputSize, int trainSize,
    bool logResgress):
m_logResgress(logResgress)
{
    m_gausProc = Vector<libgp::GaussianProcess*>(gpSize);
    m_switch = Vector<bool>(gpSize);
    m_switch.fill(false);
    m_trainCount = Vector<int>(gpSize);
    m_trainCount.fill(0);
    m_trainSize = Vector<int>(gpSize);
    m_trainSize.fill(trainSize);
    
    // Initilise all GPs
    for (int i = 0; i < gpSize; ++i)
    {
        libgp::GaussianProcess* process
            = new libgp::GaussianProcess(inputSize, "CovSEiso");
        Eigen::VectorXd params(process->covf().get_param_dim());
        for(int j = 0; j < process->covf().get_param_dim(); j++) params[j] = -1;
        process->covf().set_loghyper(params);
        m_gausProc[i] = process;
    }

    // Set deult normilisation
    Vector<double> inputNorm = Vector<double>(inputSize);
    for (int i = 0; i < inputSize; ++i) inputNorm[i] = 1.0;
    if (m_logResgress)
    {
        setNormParams(inputNorm, 10);
    } else
    {
        setNormParams(inputNorm, 1.0);
    }
    m_optimiser.init();
}

GaussianProcess::GaussianProcess(std::string gpDir)
{
    // Open meta file
    std::ifstream meta(gpDir + "/meta.dat");
    if (!meta)
    {
        // Cant find meta
        std::cerr << "Error: Cannot open meta file \"" << gpDir + "/meta"
            << "\"." << std::endl;
        exit(1);
    }
    // load meta data
    int inputSize, gpSize;
    meta >> gpSize;
    m_gausProc = Vector<libgp::GaussianProcess*>(gpSize);
    for (int i = 0; i < gpSize; ++i)
    {
        m_gausProc[i] = new libgp::GaussianProcess(("./" + gpDir + "/"
            + std::to_string(i) + ".gp").c_str());
    }
    m_switch = Vector<bool>(gpSize);
    for (int i = 0; i < gpSize; ++i) meta >> m_switch[i];
    m_trainCount = Vector<int>(gpSize);
    m_trainSize = Vector<int>(gpSize);
    for (int i = 0; i < gpSize; ++i) meta >> m_trainCount[i];
    for (int i = 0; i < gpSize; ++i) meta >> m_trainSize[i];
    meta >> inputSize;
    m_inputNorm = Vector<double>(inputSize);
    for (int i = 0; i < inputSize; ++i) meta >> m_inputNorm[i];
    meta >> m_outputNorm;
    meta >> m_logResgress;
}

GaussianProcess::~GaussianProcess()
{
    for (int i = 0; i < m_gausProc.size(); ++i)
    {
        delete m_gausProc[i];
    }
}

void GaussianProcess::save(std::string outputDir)
{
    // First save the meta data
    std::ofstream meta(outputDir + "/meta.dat");
    meta << m_gausProc.size() << "\n";
    for (int i = 0; i < m_gausProc.size(); ++i) meta << m_switch[i] << "\t";
    meta << "\n";
    for (int i = 0; i < m_gausProc.size(); ++i) meta << m_trainCount[i] << "\t";
    meta << "\n";
    for (int i = 0; i < m_gausProc.size(); ++i) meta << m_trainSize[i] << "\t";
    meta << "\n";
    meta << m_inputNorm.size() << "\n";
    for (int i = 0; i < m_inputNorm.size(); ++i) meta << m_inputNorm[i] << "\t";
    meta << "\n";
    meta << m_outputNorm << "\n";
    meta << m_logResgress;
    meta.close();

    // Save the GPs
    for (int i = 0; i < m_gausProc.size(); ++i)
    {
        m_gausProc[i]->write((outputDir + "/" + std::to_string(i)
            + ".gp").c_str());
    } 
}

void GaussianProcess::setNormParams(const Vector<double>& inputNorm,
    double outputNorm)
{
    m_inputNorm = inputNorm;
    if (m_logResgress)
    {
        m_outputNorm = std::log10(outputNorm);
    } else
    {
        m_outputNorm = outputNorm;
    }
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
        for (int i = 0; i < input.size(); ++i)
        {
            input[i] = input[i] / m_inputNorm[i];
        }
        if (m_logResgress)
        {
            output[0] = std::pow(10, m_gausProc[gpID]->f(input.begin()) * m_outputNorm);
            output[1] = m_gausProc[gpID]->var(input.begin());
        } else
        {
            output[0] = m_gausProc[gpID]->f(input.begin())
                * m_outputNorm;
            output[1] = m_gausProc[gpID]->var(input.begin());
        }
    }
}

void GaussianProcess::addData(int gpID, const Vector<double>& input,
    double output)
{
    /* Normilise data and add it to GP */
    Vector<double> inputNormed(input.size());
    for (int i = 0; i < input.size(); ++i)
    {
        inputNormed[i] = input[i] / m_inputNorm[i];
    }
    if (m_logResgress)
    {
        m_gausProc[gpID]->add_pattern(inputNormed.begin(), std::log10(output)
            / m_outputNorm);
    } else
    {
        m_gausProc[gpID]->add_pattern(inputNormed.begin(), output
            / m_outputNorm);
    }
    m_trainCount[gpID]++;
    // If training set hits the max limit then train the GP
    if (m_trainCount[gpID] == m_trainSize[gpID])
    {
        train(gpID);
    }
}

void GaussianProcess::train(int gpID)
{
    std::ofstream test("./testGp.dat");
    std::cout << "Optimising GP " << gpID << "..." << std::endl;
    m_optimiser.maximize(m_gausProc[gpID], 100, 1);
    m_trainCount[gpID] = 0; // Reset training count
    m_switch[gpID] = true;
    m_trainSize[gpID] *= 2; // twice as many points needed now
    std::cout << "Optimisation complete!" << std::endl;
}