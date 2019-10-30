#ifndef GaussianProcess_HH
#define GaussianProcess_HH

#include "gp.h"
#include "gp_utils.h"
#include "rprop.h"
#include <Eigen/Dense>


class GaussianProcess
{
public:
    /* Constuctor for a new GP */
    GaussianProcess(int gpSize, int inputSize, int trainSize);

    /* Constructor for a pre traine GP */
    GaussianProcess(std::string gpDir);

    ~GaussianProcess();

    /* Sets the normilsation partamters of the inputs and outputs */
    void setNormParams(double inputNorm[], outputNorm);

    /* Runs the GP for the given ID. returns the mean and variancve 
       in output */
    void run(int gpID, double input[], double output[2]);

    /* Adds a data point to the GP for given ID */
    void addData(int gpID, double input[], double output);

    /* Saves the GP at path outputDir */
    void save(std::string outputDir, std::string fieldName)

private:
    /* Trains the GP with given ID*/
    void train(int gpID);

private:
    libgp::GaussianProcess* m_gausProc; // Vector of GPS
    Eigen::VectorXd** m_input; // Input training data
    double** m_output;         // Output
    bool* m_switch;            // bool checking if trained
    int* m_trainCount;         // Count how many traing points saved
    int m_trainSize;           // Points aved before being trained
    int m_inputSize;             // Dimensions of input to GP
    libgp::RProp m_optimiser;  // Class to optimise GP
    double* m_inputNorm;         // Normilisation of input
    double m_outputNorm;         // Normilisation of output 
};
#endif