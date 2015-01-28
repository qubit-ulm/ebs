#ifndef LAMBDA_MAX_H
#define LAMBDA_MAX_H

#include <memory>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

namespace
{
    void initializeAvVector(const ublas::vector<double>& v,
                            ublas::vector<double>& Av)
    {
        using namespace std;
        size_t nn = v.size() - 1;

        for (size_t i = 0; i < nn; i++)
            Av[i]= v[i+1] - v[i];
    }

    /**
    * Rose algorithm to efficiently solve the
    * problem
    *
    *   A * z = b
    *
    * where
    * A=[ 2  -1   0   0;
    *	  -1   2  -1   0;
    *	   0  -1   2  -1;
    *	   0   0  -1   2]
    *
    * The function writes back it's values to vector z,
    * which is assumed to have the same size as b.
    * Additionally the maximum value of z is return in the
    * z_max value
    */
    void roseAlgorithm(const ublas::vector<double>& b,
                       ublas::vector<double>& z,
                       double & z_max)
    {
        int nn = b.size();

        double s = 0.0;
        for (int i = 0; i < nn; i++)
            s += b(i) * (i + 1);
        s /= (nn + 1);

        z(nn - 1) = b(nn - 1) - s;
        for (int i = (nn-2); i >= 0; i--) {
            z(i) = b(i) + z(i+1);
        }

        z_max = fabs(z(0));
        for (int i = 1; i < nn; i++)
        {
            z(i) += z(i-1);

            s = fabs(z(i));
            if (s > z_max)
                z_max = s;
        }
    }
}

double computeLambdaMax(const ublas::vector<double> &v)
{
    using namespace std;

    double lambda_max = -1.0;
    
    size_t n = v.size();
    ublas::vector<double>  z(n - 1);
    ublas::vector<double> Av(n - 1);

    initializeAvVector(v, Av);
    roseAlgorithm(Av, z, lambda_max);

    return lambda_max;
}

#endif // LAMBDA_MAX_H

