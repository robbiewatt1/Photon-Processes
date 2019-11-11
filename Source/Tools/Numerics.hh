#ifndef Numerics_HH
#define Numerics_HH

#include "Vector.hh"
#include "Matrix.hh"

namespace Numerics
{
    /* Gives the index of the vector closest to the value*/
    int vectorIndex(const Vector<double>& vector, double value);

    /* Simpsons method for integration */
    double simpsons(const Vector<double>& variable,
        const Vector<double>& integrand);

    /* Linear intepolation method in 1D */
    double interpolate1D(const Vector<double>& sampleX,
        const Vector<double>& sampleY, double queryX);

    /* Linear intepolation method in 2D */
    double interpolate2D(const Vector<double>& sampleX,
        const Vector<double>& sampleY, Matrix<double> sampleZ,
        double queryPoint[2]);
}
#endif