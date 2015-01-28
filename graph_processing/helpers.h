#ifndef GRAPHHELPERS_H
#define GRAPHHELPERS_H

#include <iostream>
#include <boost/graph/adjacency_list.hpp>

#include "../common/helpers.h"

// Boost Graph
namespace helpers
{
    template <class CapacityMapType,class ReverseMapType, class ResidualCapacityType>
    class EdgeWriter {

    public:
        EdgeWriter(CapacityMapType c, ReverseMapType r, ResidualCapacityType rc)
            : m_capacity_prop(c),
              m_rev_prop(r),
              m_residual_capacity_prop(rc)
        { }

        template <class EdgeDescriptorT>
        void operator()(std::ostream &out, const EdgeDescriptorT& edge_desc) const
        {
            /*
            auto capacity = get(m_capacity_prop, edge_desc);
            auto reverse = get(m_rev_prop, edge_desc);
            auto residual = get(m_residual_capacity_prop, edge_desc);
            */
            //out << "[label=\"" << get(m_capacity_prop, e) << "\"]";//, taillabel=\"" << m_rev_prop[e] << "\"]";
            out << "[label=\"" << "\"]";
        }

    private:
        CapacityMapType m_capacity_prop;
        ReverseMapType m_rev_prop;
        ResidualCapacityType m_residual_capacity_prop;
    };

    template <class CapacityMapType,class ReverseMapType, class ResidualCapacityType>
    inline EdgeWriter<CapacityMapType,ReverseMapType, ResidualCapacityType>
    make_edge_writer(CapacityMapType c, ReverseMapType r, ResidualCapacityType rc)
    {
        return EdgeWriter<CapacityMapType,ReverseMapType,ResidualCapacityType>(c, r, rc);
    }
}

// Misc Helpers
namespace helpers
{
    inline std::pair<double, double> findMinMax(const std::vector<double>& input)
    {
        auto min_max_it = std::minmax(input.begin(), input.end());
        return std::make_pair(*min_max_it.first, *min_max_it.second);
    }

    inline std::tuple<std::vector<int>, std::vector<int>> findMinMaxValueCombinationIdcs(const std::vector<double>& data,
                                                                                         const std::vector<double>& weights,
                                                                                         const std::vector<double>& labels)
    {
        using namespace std;

        vector<int> data_idcs;
        data_idcs.push_back(distance(data.begin(), // min data idx
                                     min_element(data.begin(), data.end())) + 2);
        data_idcs.push_back(distance(data.begin(), // max data idx
                                     max_element(data.begin(), data.end())) + 2);
        data_idcs.push_back(distance(weights.begin(), // max weight idx
                                     max_element(weights.begin(), weights.end())) + 2);

        vector<int> label_idcs;
        label_idcs.push_back(distance(labels.begin(), // min label idx
                                      min_element(labels.begin(), labels.end())));
        label_idcs.push_back(distance(labels.begin(), // max label idx
                                      max_element(labels.begin(), labels.end())));

        return make_tuple(data_idcs, label_idcs);
    }

    inline long long maxDataCosts(std::function<long long(std::tuple<int, int>, int)> data_cost_fn,
                                  const std::vector<std::pair<double, double>>& csv_data,
                                  const std::vector<double>& labels)
    {
        return 0;
    }

    inline long long maxDataCosts(std::function<long long(std::tuple<int, int>, int)> data_cost_fn,
                                  const std::vector<double>& data,
                                  const std::vector<double>& weights,
                                  const std::vector<double>& labels)
    {
        using namespace std;

        vector<int> data_idcs;
        vector<int> label_idcs;
        tie(data_idcs, label_idcs) = findMinMaxValueCombinationIdcs(data, weights, labels);

        vector<long long> costs;
        for (auto data_idx : data_idcs) {
            for (auto label_idx : label_idcs) {
                costs.push_back(data_cost_fn(make_tuple(data_idx, label_idx), 0));
            }
        }

        auto max_cost_it = max_element(costs.begin(), costs.end());
        return *max_cost_it;
    }

    inline long long maxSmoothCosts(std::function<long long(std::tuple<int, int, int, int>, int)> smooth_cost_fn,
                                    const std::vector<std::pair<double, double>>& csv_data,
                                    const std::vector<double>& labels)
    {
        return 0;
    }

    inline long long maxSmoothCosts(std::function<long long(std::tuple<int, int, int, int>, int)> smooth_cost_fn,
                                    const std::vector<double>& data,
                                    const std::vector<double>& weights,
                                    const std::vector<double>& labels)
    {
        using namespace std;

        vector<int> data_idcs;
        vector<int> label_idcs;
        tie(data_idcs, label_idcs) = findMinMaxValueCombinationIdcs(data, weights, labels);

        vector<long long> costs;

        for (size_t i = 0; i < labels.size() / 2; i++) {
            for (size_t j = labels.size() / 2; j < labels.size(); j++) {
                costs.push_back(smooth_cost_fn(make_tuple(data_idcs[2], data_idcs[2], i, j), 0));
            }
        }

        auto max_cost_it = max_element(costs.begin(), costs.end());
        return *max_cost_it;
    }

    inline long long discretizeAndReweightCost(double cost, long long max_cost)
    {
        if (max_cost <= 0)
            return (long long)cost;

        return (long long)((cost * 10000.0) / (double)max_cost) ;
    }

    inline long long normalizeToOneAndDiscretize(double cost)
    {
        using namespace std;

        //double y = 1.0 - 1.0 / (1.0 + cost / 1000.0);
        //return (long long)(y * 1000.0);

        //double y = tanh(cost);
        //return (long long)(y * 1000.0);

        double y = 1.0 / (1.0 + exp(-1.0 * cost));
        return (long long)((2 * y - 1) * 1000.0);
    }
}

#endif // GRAPHHELPERS_H
