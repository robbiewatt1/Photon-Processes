#ifndef Numerics_HH
#define Numerics_HH

namespace Numerics
{

    /* Basic simpsons method for integration */
    double simpsons(double* variable, double* integrand, int resolusion);

    /* Basic linear intepolation method in 1D */
    double interpolate1D(double* sampleX, double* sampleY,
            int sampleSize, double queryX);

    /* Basic linear intepolation method in 2D */
    double interpolate2D(double* sampleX, double* sampleY, double** sampleZ,
            int sampleSize[2], double queryPoint[2]);
}
#endif