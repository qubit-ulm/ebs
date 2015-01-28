#ifndef TUPLE_HELPER_H
#define TUPLE_HELPER_H

#include <iostream>
#include <cmath>

namespace helpers
{
    template <typename VectorType>
    void diff(const VectorType &v,
    	      VectorType &dv)
    {
        using namespace std;

        size_t n = v.size();
        dv.resize(n - 1);

        for (size_t i = 1; i < n; i++)
            dv(i-1) = v(i) - v(i-1);
    }

    template <typename VectorType>
    double countJumpsInDiff(const VectorType &dv,
                            double thresh)
    {
        using namespace std;

        size_t n = dv.size();
        double s = 1.0;

        for (size_t i = 0; i < n; i++) {
            auto dv_i = fabs(dv(i));
            s += dv_i > thresh ? 1.0 : 0.0;
        }

        return s;
    }

    template <typename VectorType>
    void populateDataAndWeightVectors(const VectorType &input,
                                      const VectorType &dv,	
                                      VectorType &data,
                                      VectorType &weights,
			    	                  double thresh = 0.0)
    {
        using namespace std;

        size_t  n = dv.size();
        size_t nn = countJumpsInDiff(dv, thresh);
        data.resize(nn);
        weights.resize(nn);

        data(0) = input(0);

        size_t j = 1;
        size_t prev_i = 0;
        for (size_t i = 1; i < n; i++) 
        {
            if (fabs(dv(i)) <= thresh)
                continue;

            data(j) = input(i);
            weights(j-1) = i - prev_i;
            j++;

            prev_i = i;
        }

        weights(nn- 1) = n - prev_i + 1;
    }

    template <typename VectorType>
    void combineConsecutiveValues(VectorType &data, 
                                  VectorType &weights)
    {
        using namespace std;
        
        VectorType cp_data(data);
        VectorType cp_weights(weights);

        size_t n = data.size();
        size_t j = 0;
        double prev_data = sqrt(-1.0);
        
        for (size_t i = 0; i < n; i++)
        {
            double cur_data = cp_data(i);

            if (! std::isnan(prev_data) && prev_data != cur_data) 
            {
                j = j + 1;
                data(j) = cur_data;
                weights(i) = 0;
            }

            data(j) = cur_data;
            weights(j) = weights(j) + cp_weights(i);

            prev_data = cur_data;
        }
    }

    template <typename VectorType>
    void postprocessTVDNData(const VectorType &input,
		    	             VectorType &data,
			                 VectorType &weights)
    {
	    using namespace std;

	    VectorType dv;
	    diff(input, dv);
	    populateDataAndWeightVectors(input, dv, data, weights);
        //combineConsecutiveValues(data, weights);
    }
}

#endif // TUPLE_HELPER_H
