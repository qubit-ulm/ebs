#ifndef LAMBDA_OPT_H
#define LAMBDA_OPT_H

#include <iomanip>
#include <memory>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/log/trivial.hpp>

#include "../common/tuple_helper.h"
#include "condat_denoise.h"

namespace ublas = boost::numeric::ublas;

namespace {

    struct BisectionMethod
    {
        typedef ublas::vector<double> Vec;

        BisectionMethod(const Vec &v)
            : m_noisy(v)
            , m_denoised(v.size())
            , m_diff(v.size() - 1)
        { }

        double findLambdaOpt(double lambda_min,
                             double lambda_max,
                             double n_min = -1,
                             double n_max = -1,
                             std::size_t remaining_iters = 50)
        {
            using namespace std;

            double thr = 1e-7;
            double lambda_pivot = (lambda_min + lambda_max) / 2;

            if (remaining_iters == 0)
	    {
		BOOST_LOG_TRIVIAL(debug) << "Max iterations reached, "
					    "returning lambda_pivot: "
					 << lambda_pivot;
                return lambda_pivot;
	    }

            n_min = n_min < 0 ? countJumpsForLambda(lambda_min, thr) : n_min;
            n_max = n_max < 0 ? countJumpsForLambda(lambda_max, thr) : n_max;
            double n_pivot = countJumpsForLambda(lambda_pivot, thr);

            double s_1 = slope(n_min, n_pivot, lambda_min, lambda_pivot);
            double s_2 = slope(n_pivot, n_max, lambda_pivot, lambda_max);

	    BOOST_LOG_TRIVIAL(debug) 
		<< scientific 
		<< setw(12) << "lambda_min: " << lambda_min << " "
		<< setw(12) << "lambda_pvt: " << lambda_pivot << " "
		<< setw(12) << "lambda_max: " << lambda_max << " "
		<< setw(8) << "n_min: " << n_min << " "
		<< setw(8) << "n_pvt: " << n_pivot << " "
		<< setw(8) << "n_max: " << n_max << " "
		<< setw(8) << "s_1: " << s_1 << " "
		<< setw(8) << "s_2: " << s_2;

            if (s_1 > s_2)
                return findLambdaOpt(lambda_min, lambda_pivot, 
				     n_min, n_pivot, remaining_iters - 1);
            else
                return findLambdaOpt(lambda_pivot, lambda_max, 
				     n_pivot, n_max, remaining_iters - 1);
        }

    private:
        const Vec &m_noisy;
        Vec m_denoised;
        Vec m_diff;

    private:
        void diff(const ublas::vector<double> &v)
        {
	    	helpers::diff(v, m_diff);
        }

        double countJumpsInDiff(double thresh)
        {
	   		return helpers::countJumpsInDiff(m_diff, thresh);
		}

        double slope(const double n_1,
                     const double n_2,
                     const double l_1,
                     const double l_2)
        {
            double s = (n_2 - n_1) / (l_2 - l_1);
            return fabs(s);
        }

        double countJumpsForLambda(const double lambda,
                                   const double thresh)
        {
            TV1D_denoise(m_noisy, m_denoised, lambda);
            diff(m_denoised);

            auto n = countJumpsInDiff(thresh);
		    return n;
        }
    };

    struct SteepDecentMethod
    {
        typedef ublas::vector<double> Vec;

        SteepDecentMethod(const Vec &v)
            : m_noisy(v)
            , m_denoised(v.size())
            , m_diff(v.size() - 1)
        { }

		double findLambdaOpt(double lambda_max,
							 std::size_t max_iter = 50)
		{
			using namespace std;

			auto thresh = 1.0e-7;
			auto rho = 5.0;
			auto N = (double)m_noisy.size();
			auto f_prev = 1.0;
			auto n_prev = countJumpsForLambda(f_prev * lambda_max, thresh);
			auto f = f_prev / 2.0;
			auto n = countJumpsForLambda(f * lambda_max, thresh);

			auto start_slope = calculateSlope(0.0, 1.0, N, n_prev,
											  lambda_max, thresh);

			BOOST_LOG_TRIVIAL(debug) << scientific
									 << "start_slope: "
									 << start_slope;
			
			for (size_t i = 0; i < max_iter; i++)
			{
				auto slope = calculateSlope(f, f_prev, n, n_prev,
											lambda_max, thresh); 

				BOOST_LOG_TRIVIAL(debug) 
					<< setw(12) << "f: " << f << " "
					<< setw(12) << "slope: " << slope;

				if (slope > start_slope)
					break;
				
				f_prev = f;
				n_prev = n;
				f = f_prev / rho;
				n = countJumpsForLambda(f * lambda_max, thresh);
			}

			double lambda_opt = f * lambda_max;
			return lambda_opt;
		}

    private:
		const Vec &m_noisy;
		Vec m_denoised;
		Vec m_diff;

    private:
		void diff(const ublas::vector<double> &v)
        {
	    	helpers::diff(v, m_diff);
        }

        double countJumpsInDiff(double thresh)
        {
	   		return helpers::countJumpsInDiff(m_diff, thresh);
		}

		double countJumpsForLambda(const double lambda,
								   const double thresh)
		{
			using namespace std;

			TV1D_denoise(m_noisy, m_denoised, lambda);
			diff(m_denoised);

			auto N = (double)m_noisy.size();
			auto n = countJumpsInDiff(thresh);
			return min(n, N);
		}

		double calculateSlope(double f_min, 
							  double f_max,
							  double n_min,
							  double n_max,
							  double lambda_max,
							  const double thresh)
		{
			auto slope = fabs(n_min - n_max) / (f_max - f_min);
			return slope;
		}
	};
}

double computeLambdaOpt(const BisectionMethod::Vec &v,
                        double lambda_min,
                        double lambda_max)
{
    using namespace std;

    //BisectionMethod algo(v);
    //double lambda_opt = algo.findLambdaOpt(lambda_min, lambda_max);
    
	SteepDecentMethod algo(v);
	double lambda_opt = algo.findLambdaOpt(lambda_max);
	
	return lambda_opt;
}

#endif // LAMBDA_OPT_H
