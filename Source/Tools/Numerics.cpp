#include "Numerics.hh"
#include <iostream>

int Numerics::vectorIndex(const Vector<double>& vector, double value)
{
    int index;
    if (value < vector[0])
    {
        index = 0;
    } else if (value > *vector.end())
    {
        index = vector.size() - 1;
    } else
    {
        for (int i = 0; i < vector.size(); i++)
        {
            if (vector[i] > value)
            {
                double closest = vector[i] + vector[i-1]
                        - 2.0 * value;
                index = closest > 0 ? (i - 1) : i;
                break;
            }
        }
    }
    return index;
}

double Numerics::simpsons(const Vector<double>& variable,
    const Vector<double>& integrand)
{
    double deltaX = (*variable.end() - variable[0]) / (variable.size());
    double integral = integrand[0];
    for (int i = 1; i < variable.size() - 1; i += 2)
    {
        integral += 4.0 * integrand[i];
    }
    for (int i = 2; i < variable.size() - 2; i += 2)
    {
        integral += 2.0 * integrand[i];
    }
    integral += *integrand.end();

    return (integral * deltaX / 3.0);
}

double Numerics::interpolate1D(const Vector<double>& sampleX,
    const Vector<double>& sampleY, double queryX)
{
    // Get the closest points
    int sampleSize = sampleX.size();
    int closestIndex[2];
    if (queryX < sampleX[0])
    {
        std::cout << "Warning: Extrapolation used in interpolate2D()"
            << std::endl;
        closestIndex[0] = 0;
        closestIndex[1] = 1;
    } else if (queryX > *sampleX.end())
    {
        std::cout << "Warning: Extrapolation used in interpolate2D()"
            << std::endl;
        closestIndex[0] = sampleSize - 2;
        closestIndex[1] = sampleSize - 1;      
    } else
    {
        for (int i = 0; i < sampleSize; i++)
        {
            if(sampleX[i] > queryX)
            {
                closestIndex[0] = i - 1;
                closestIndex[1] = i;
                break;
            }
        }
    }
    double x1 = sampleX[closestIndex[0]];
    double x2 = sampleX[closestIndex[1]];
    double y1 = sampleY[closestIndex[0]];
    double y2 = sampleY[closestIndex[1]];
    return 1.0 / (x2 - x1) * (y1 * (x2 - queryX) + y2 * (queryX - x1));
}

double Numerics::interpolate2D(const Vector<double>& sampleX,
    const Vector<double>& sampleY, Matrix<double> sampleZ,
    double queryPoint[2])
{
    // Get the closest points
    int sampleSize[2] = {sampleX.size(), sampleY.size()};
    int closestIndexX[2];
    if (queryPoint[0] < sampleX[0])
    {
        std::cout << "Warning: Extrapolation used in interpolate2D()"
            << std::endl;
        closestIndexX[0] = 0;
        closestIndexX[1] = 1;
    } else if (queryPoint[0] > *sampleX.end())
    {
        std::cout << "Warning: Extrapolation used in interpolate2D()"
            << std::endl;
        closestIndexX[0] = sampleSize[0] - 2;
        closestIndexX[1] = sampleSize[0] - 1;      
    } else
    {
        for (int i = 0; i < sampleSize[0]; i++)
        {
            if(sampleX[i] > queryPoint[0])
            {
                closestIndexX[0] = i - 1;
                closestIndexX[1] = i;
                break;
            }
        }
    }
    int closestIndexY[2];
    if (queryPoint[1] < sampleY[0])
    {
        std::cout << "Warning: Extrapolation used in interpolate2D()"
            << std::endl;
        closestIndexY[0] = 0;
        closestIndexY[1] = 1;
    } else if (queryPoint[1] > *sampleY.end())
    {
        std::cout << "Warning: Extrapolation used in interpolate2D()"
            << std::endl;
        closestIndexY[0] = sampleSize[1] - 2;
        closestIndexY[1] = sampleSize[1] - 1;
    } else
    {
        for (int i = 0; i < sampleSize[1]; i++)
        {
            if(sampleY[i] > queryPoint[1])
            {
                closestIndexY[0] = i - 1;
                closestIndexY[1] = i;
                break;
            }
        }
    }
    double x1 = sampleX[closestIndexX[0]];
    double x2 = sampleX[closestIndexX[1]];
    double y1 = sampleY[closestIndexY[0]];
    double y2 = sampleY[closestIndexY[1]];
    double z11 = sampleZ[closestIndexX[0]][closestIndexY[0]];
    double z12 = sampleZ[closestIndexX[0]][closestIndexY[1]];
    double z21 = sampleZ[closestIndexX[1]][closestIndexY[0]];
    double z22 = sampleZ[closestIndexX[1]][closestIndexY[1]];  
    return 1.0 / ((x2 - x1) * (y2 - y1))
        * (z11 * (x2 - queryPoint[0]) * (y2 - queryPoint[1])
         + z21 * (queryPoint[0] - x1) * (y2 - queryPoint[1])
         + z12 * (x2 - queryPoint[0]) * (queryPoint[1] - y1)
         + z22 * (queryPoint[0] - x1) * (queryPoint[1] - y1));
}