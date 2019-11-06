#ifndef Numerics_HH
#define Numerics_HH

#include "Vector.hh"
#include "Matrix.hh"

namespace Numerics
{
    /* Basic simpsons method for integration */
    double simpsons(const Vector<double>& variable,
        const Vector<double>& integrand);

    /* Basic linear intepolation method in 1D */
    double interpolate1D(const Vector<double>& sampleX,
        const Vector<double>& sampleY, double queryX);

    /* Basic linear intepolation method in 2D */
    double interpolate2D(const Vector<double>& sampleX,
        const Vector<double>& sampleY, Matrix<double> sampleZ,
        double queryPoint[2]);
}
#endif