#ifndef BINOPT_H
#define BINOPT_H

#include <functional>
#include <map>
#include <iomanip>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

#include <boost/log/trivial.hpp>

#include "energy.h"
#include "sitesstore.h"
#include "runtime_statistics.h"

using namespace std;

template<typename EnergyType = long long,
         typename DataCostFnArgType = int,
         typename SmoothCostFnArgType = int,
         typename LabelCostFnArgType = int>
class BinaryOptimization
{
public:
    typedef typename EnergyGraph<EnergyType>::VertexDescriptor VertexDescriptor;
    typedef typename EnergyGraph<EnergyType>::EdgeDescriptor EdgeDescriptor;

public:

    BinaryOptimization(int n_sites, int n_labels)
        : m_last_expansion_energy(numeric_limits<EnergyType>::max())
        , m_record_energy_graph_dumps(false)
        , m_record_energy_history(true)
    {
        vector<VertexDescriptor> vertex_descs;

        initializeEnergyGraph(n_sites, vertex_descs);
        initializeSitesStore(vertex_descs);
        initializeLabelTable(n_labels);
    }

    ~BinaryOptimization()
    { }

    EnergyType initiallyAssignLabels()
    {
        vector<VertexDescriptor> vertices;
        m_sites_store.queryAllVertices(vertices);

        //auto initial_cost = initiallyAssignLabelsByMinDataCost(vertices);
        auto initial_cost = initallyAssignLabelsRandomly(vertices);

        updateCountingStatistics();
        return initial_cost;
    }

    EnergyType expansion(int max_iterations = 100)
    {
        BOOST_LOG_TRIVIAL(debug) << "*******************************************";
        BOOST_LOG_TRIVIAL(debug) << "* starting alpha expansion";
        BOOST_LOG_TRIVIAL(debug) << "*";

        EnergyType new_energy = -1;

        if (max_iterations < 0) {
            new_energy = expansionConcentratingOnEnergyReducingLabels();
        }
        else {
            new_energy = expansionSweepingAllLabels(max_iterations);
        }

        BOOST_LOG_TRIVIAL(debug) << "*";
        BOOST_LOG_TRIVIAL(debug) << "*******************************************";
        BOOST_LOG_TRIVIAL(debug) << "";

        return new_energy;
    }

    BinaryOptimization& setDataCost(function<EnergyType(tuple<int, int>, DataCostFnArgType)> data_cost_fn)
    {
        m_data_cost_fn = data_cost_fn;
        return *this;
    }

    BinaryOptimization& setSmoothnessCost(function<EnergyType(tuple<int, int, int, int>, SmoothCostFnArgType)> smooth_cost_fn)
    {
        m_smooth_cost_fn = smooth_cost_fn;
        return *this;
    }

    BinaryOptimization& setLabelCost(function<EnergyType (tuple<int, int, int, int>, LabelCostFnArgType)> label_cost_fn)
    {
        m_label_cost_fn = label_cost_fn;
        return *this;
    }

    int whichLabel(VertexDescriptor site_id)
    {
        return m_sites_store.whichLabel(site_id);
    }

    vector<int> whichLabels()
    {
        vector<VertexDescriptor> vertices;
        m_sites_store.queryAllVertices(vertices);

        vector<int> assignments(vertices.size());

        for (auto v : vertices)
        {
            auto idx = whichVertexIndex(v);
            assignments[idx] = m_sites_store.whichLabel(v);
        }

        return assignments;
    }

    void recordEnergyGraphDumps(bool record_dumps = true)
    {
        m_record_energy_graph_dumps = record_dumps;
    }

    void recordEnergyHistory(bool record_history = true)
    {
        m_record_energy_history = record_history;
    }

    std::vector<EnergyType>& findEnergyHistory(std::string label)
    {
        return m_runtime_statistics.findEnergyHistory(label);
    }

private:
    vector<int> m_label_table;
    //vector<int> m_label_costs;
    EnergyGraph<EnergyType> m_energy_graph;

    function<EnergyType (tuple<int, int>, DataCostFnArgType)> m_data_cost_fn;
    function<EnergyType (tuple<int, int, int, int>, SmoothCostFnArgType)> m_smooth_cost_fn;
    function<EnergyType (tuple<int, int, int, int>, LabelCostFnArgType)> m_label_cost_fn;

    SitesStore<EnergyType, VertexDescriptor, int> m_sites_store;

    EnergyType m_last_expansion_energy;
    RuntimeStatistics<std::string, EnergyType> m_runtime_statistics;

    bool m_record_energy_graph_dumps;
    bool m_record_energy_history;

private:

    void initializeEnergyGraph(const int n_sites,
                               vector<VertexDescriptor>& vertex_descs)
    {
        vertex_descs.clear();

        for (int i = 0; i < n_sites; i++)
        {
            auto v = m_energy_graph.addVariable();
            if (i > 0) {
                m_energy_graph.addTerm2(v, vertex_descs.back(),
                                        0, 0, 0, 0);
            }

            vertex_descs.push_back(v);
        }
    }

    EnergyType initallyAssignLabelsRandomly(vector<VertexDescriptor>& vertices)
    {
        if (! m_data_cost_fn)
            return -1;

        permuteLabelTable();

        EnergyType initial_cost = 0.0;
        int n = 0;
        for (auto vertex_desc : vertices)
        {
            auto vertex_idx = whichVertexIndex(vertex_desc);
            auto label = m_label_table[n % m_label_table.size()];
            n++;
            auto cost = m_data_cost_fn(make_tuple(vertex_idx, label), 0);

            m_sites_store.assignLabel(vertex_desc, label);
            m_sites_store.assignDataCost(vertex_desc, cost);

            initial_cost += cost;
        }

        return -1;
    }

    EnergyType initiallyAssignLabelsByMinDataCost(vector<VertexDescriptor>& vertices)
    {
        if (! m_data_cost_fn)
            return -1;

        EnergyType initial_cost = 0.0;
        for (auto vertex_desc : vertices)
        {
            auto best_label = findMinDataCostLabel(vertex_desc);
            auto cost = m_data_cost_fn(make_tuple(vertex_desc, best_label), 0);

            m_sites_store.assignLabel(vertex_desc, best_label);
            m_sites_store.assignDataCost(vertex_desc, cost);

            initial_cost += cost;
        }

        return initial_cost;
    }

    void initializeSitesStore(vector<VertexDescriptor>& vertex_descs)
    {
        for (auto v : vertex_descs)
        {
            m_sites_store.addVertex(v)
                         .setActive(v, true);
        }
    }

    void initializeLabelTable(int n_labels)
    {
        for (int i = 0; i < n_labels; i++)
            m_label_table.push_back(i);
    }

    int findMinDataCostLabel(VertexDescriptor vertex_desc)
    {
        auto min_cost = std::numeric_limits<EnergyType>::max();
        int min_label = 0;

        for (int cur_label : m_label_table) {
            auto cost = m_data_cost_fn(make_tuple(vertex_desc, cur_label), 0);
            if (cost >= min_cost)
                continue;

            min_cost = cost;
            min_label = cur_label;
        }

        return min_label;
    }

    void proposeAlphaLabel(int& next_label,
                           int& cycle_size,
                           int& cycle)
    {
        auto alpha_label = m_label_table[next_label];
        if (alphaExpansion(cycle, 0, alpha_label)) {
            // keep label for smaller queue, because it lead to energy change
            next_label++;
            updateLabelInformation();
        }
        else {
            // don't keep the label, because it is uninteressting
            cycle_size--;
            swap(m_label_table[next_label], m_label_table[cycle_size]);
        }
    }

    void adjustLabelQueue(queue<int>& sizes_queue,
                          int& next_label,
                          int& start_label,
                          int& cycle_size)
    {
        if (next_label == start_label) {
            // No expansion was successfull, to try more from the previous queue
            next_label = sizes_queue.back();
            sizes_queue.pop();
        }
        else if (cycle_size < (int)(sizes_queue.back() / 2.0)) {
            // Some expansion were successfull, so focus on them
            next_label = 0;
            sizes_queue.push(cycle_size);
        }
        else {
            // All expansion where successfull, so complete sweep
            next_label = 0;
        }
    }

    EnergyType expansionConcentratingOnEnergyReducingLabels()
    {
        queue<int> sizes_queue;
        sizes_queue.push(m_label_table.size());

        permuteLabelTable();
        updateLabelInformation();

        int next_label = 0;

        for (auto cycle = 0; !sizes_queue.empty(); cycle++)
        {
            // Pass over unchecked labels in the current queue
            int start_label = next_label;
            int cycle_size = sizes_queue.back();
            int steps_in_cycle = cycle_size - start_label;

            do {
                proposeAlphaLabel(next_label, cycle_size, cycle);
            } while (next_label < cycle_size);

            adjustLabelQueue(sizes_queue, next_label, start_label, cycle_size);
        }

        return computeEnergy();
    }

    EnergyType expansionSweepingAllLabels(int max_iterations)
    {
        EnergyType new_energy = computeEnergy();
        EnergyType old_energy = new_energy + 1;

        for (int i = 0; i < max_iterations; i++)
        {
            old_energy = new_energy;
            new_energy = doExpansionIteration(i);

            if (new_energy == old_energy)
                break;
        }

        return new_energy;
    }

    EnergyType doExpansionIteration(int iter)
    {
        updateLabelInformation();
        permuteLabelTable();

        int label_iter = 0;
        for (auto label : m_label_table)
         {
            BOOST_LOG_TRIVIAL(debug) << "\t----------------------------";
            BOOST_LOG_TRIVIAL(debug) << "\tIter: " << iter;
            BOOST_LOG_TRIVIAL(debug) << "\tAttempting label: " << label;

            alphaExpansion(iter, label_iter, label);
            label_iter++;
        }

        return computeEnergy(); // compute energy
    }

    void updateLabelInformation()
    {
        //updateCountingStatistics();
        //updateLabelActiveFlags();
        updateDataCostsBasedOnCurrentLabeling();
    }

    inline void updateCountingStatistics()
    {
        m_sites_store.updateCountingStatistics();
    }

    inline void updateLabelActiveFlags()
    {
        m_sites_store.markAllVerticesInactive();
        for (auto l : m_label_table) {
            auto count = m_sites_store.labelCount(l);
            if (count <= 0)
                continue;

            m_sites_store.setActiveForLabel(l, true);
        }
    }

    inline void updateDataCostsBasedOnCurrentLabeling()
    {
        if (! m_data_cost_fn)
            return;

        vector<VertexDescriptor> vertices;
        m_sites_store.queryAllVertices(vertices);

        for (auto vertex_desc : vertices)
        {
            auto cur_label = m_sites_store.whichLabel(vertex_desc);
            auto vertex_idx = whichVertexIndex(vertex_desc);
            auto cost = m_data_cost_fn(make_tuple(vertex_idx, cur_label), 0);

            m_sites_store.assignDataCost(vertex_desc, cost);
        }
    }

    void permuteLabelTable()
    {
        srand(time(0));
        random_shuffle(m_label_table.begin(), m_label_table.end());
    }

    bool alphaExpansion(int iter, int label_iter, int alpha_label)
    {
        // Get list of active sites based on the alpha_label
        vector<VertexDescriptor> active_sites;
        m_sites_store.queryAllVertices(active_sites);
        if (active_sites.size() == 0)
        {
            BOOST_LOG_TRIVIAL(debug) << "\tNo actives vertices, skipping alpha expansion";
            return false;
        }

        // Create binary variables for each remaining site, add data costs
        // and compute the smooth costs between variables
        m_energy_graph.recycle();
        addDataCostEdges(alpha_label, active_sites, m_energy_graph);
        addSmoothingCostEdges(alpha_label, active_sites, m_energy_graph);
        addLabelCostEdges(alpha_label, active_sites, m_energy_graph);

        EnergyType energy_after_expansion = m_energy_graph.minimize();
        BOOST_ASSERT(energy_after_expansion >= 0);

        BOOST_LOG_TRIVIAL(debug) << "Energy after expansion: " << energy_after_expansion
                                 << ",\tprev expansion Energy: " << m_last_expansion_energy;

        dumpEnergyGraph(iter, label_iter, alpha_label, energy_after_expansion);
        recordEnergyHistory(iter, label_iter, alpha_label, energy_after_expansion);

        bool is_energy_improved = energy_after_expansion < m_last_expansion_energy;
        if (is_energy_improved)
        {
            acceptNewLabeling(m_energy_graph, alpha_label, active_sites);
            updateLabelInformation();
            m_last_expansion_energy = energy_after_expansion;
        }

        return is_energy_improved;
    }

    void addDataCostEdges(const int alpha_label,
                          const vector<VertexDescriptor>& active_vertices,
                          EnergyGraph<EnergyType>& energy)
    {
        if (! m_data_cost_fn)
            return;

        for (auto vertex_desc : active_vertices)
        {
            auto vertex_idx = whichVertexIndex(vertex_desc);

            // TODO: Add correct handling of function arg here
            auto e0 = m_sites_store.dataCost(vertex_desc);
            auto args = make_tuple(vertex_idx, alpha_label);
            auto e1 = safeInvokeCostFn(m_data_cost_fn, args, 0);

            energy.addTerm1(vertex_desc, e0, e1);
        }
    }

    void addSmoothingCostEdges(const int alpha_label,
                               const vector<VertexDescriptor>& active_vertices,
                               EnergyGraph<EnergyType>& energy)
    {
        if (! m_smooth_cost_fn)
            return;

        addSmoothingTypeCostEdges(m_smooth_cost_fn,
                                  alpha_label,
                                  active_vertices,
                                  energy);
    }

    void addLabelCostEdges(const int alpha_label,
                           const vector<VertexDescriptor>& active_vertices,
                           EnergyGraph<EnergyType>& energy)
    {
        using namespace std;

        if (! m_label_cost_fn)
            return;

        addSmoothingTypeCostEdges(m_label_cost_fn,
                                  alpha_label,
                                  active_vertices,
                                  energy);
    }

    inline bool isActiveNeighbour(const vector<VertexDescriptor>& actives,
                                  VertexDescriptor neighb)
    {
        return find(actives.begin(), actives.end(), neighb) != actives.end();
    }


    template<typename FnType>
    void addSmoothingTypeCostEdges(FnType cost_fn,
                                   const int alpha_label,
                                   const vector<VertexDescriptor>& active_vertices,
                                   EnergyGraph<EnergyType>& energy)
    {
        if (! cost_fn)
            return;

        for (auto vert_desc : active_vertices)
        {
            auto neighbouring_vertices = energy.neighboursOf(vert_desc);
            for (auto nb_vert_desc : neighbouring_vertices)
            {
                if (isActiveNeighbour(active_vertices, nb_vert_desc))
                {
                    addSmoothingTypeCostsForActiveNeighbourEdge(cost_fn,
                                                                alpha_label,
                                                                vert_desc,
                                                                nb_vert_desc,
                                                                energy);
                }
                else
                {
                    addSmoothingTypeCostsForInactiveNeighbourEdge(cost_fn,
                                                                  alpha_label,
                                                                  vert_desc,
                                                                  nb_vert_desc,
                                                                  energy);
                }
            }
        }
    }

    template<typename FnType>
    inline void addSmoothingTypeCostsForActiveNeighbourEdge(const FnType cost_fn,
                                                            const int alpha_label,
                                                            const VertexDescriptor& vert_desc,
                                                            const VertexDescriptor& nb_vert_desc,
                                                            EnergyGraph<EnergyType>& energy)
    {
        auto cur_label = m_sites_store.whichLabel(vert_desc);
        auto nb_label = m_sites_store.whichLabel(nb_vert_desc);

        auto vert_idx = whichVertexIndex(vert_desc);
        auto nb_idx = whichVertexIndex(nb_vert_desc);

        // TODO: Add correct handling of function arg here
        auto args = make_tuple(vert_idx, nb_idx, alpha_label, alpha_label);
        auto e00 = safeInvokeCostFn(cost_fn, args, 0);

        args = make_tuple(vert_idx, nb_idx, alpha_label, nb_label);
        auto e01 = safeInvokeCostFn(cost_fn, args, 0);

        args = make_tuple(vert_idx, nb_idx, cur_label, alpha_label);
        auto e10 = safeInvokeCostFn(cost_fn, args, 0);

        args = make_tuple(vert_idx, nb_idx, cur_label, nb_label);
        auto e11 = safeInvokeCostFn(cost_fn, args, 0);

        if (e00 + e11 > e01 + e10)
        {
            //BOOST_LOG_TRIVIAL(warning) << "Non submoduler smoothness term detected";
            //return;
            healSubmodularEnergies(e00, e01, e10, e11);
        }

        energy.addTerm2(vert_desc, nb_vert_desc, e00, e01, e10, e11);
    }

    inline void healSubmodularEnergies(EnergyType& e00,
                                       EnergyType& e01,
                                       EnergyType& e10,
                                       EnergyType& e11)
    {

        for (int i = 0; e00 + e11 > e01 + e10; i++)
        {
            if (i % 3 == 0)
                e01++;
            else if (i % 3 == 1)
                e10++;
            else
                e00--;
        }
    }

    template<typename FnType>
    inline void addSmoothingTypeCostsForInactiveNeighbourEdge(const FnType cost_fn,
                                                              const int alpha_label,
                                                              const VertexDescriptor& vert_desc,
                                                              const VertexDescriptor& nb_vert_desc,
                                                              EnergyGraph<EnergyType>& energy)
    {
        auto cur_label = m_sites_store.whichLabel(vert_desc);
        auto nb_label = m_sites_store.whichLabel(nb_vert_desc);

        auto vert_idx = whichVertexIndex(vert_desc);
        auto nb_idx = whichVertexIndex(nb_vert_desc);

        auto args = make_tuple(vert_idx, nb_idx, alpha_label, nb_label);
        auto e0 = safeInvokeCostFn(cost_fn, args, 0);

        args = make_tuple(vert_idx, nb_idx, cur_label, nb_label);
        auto e1 = safeInvokeCostFn(cost_fn, args, 0);

        energy.addTerm1(vert_desc, e0, e1);
    }

    void acceptNewLabeling(EnergyGraph<EnergyType>& energy,
                           const int alpha_label,
                           const vector<VertexDescriptor>& active_sites)
    {
        if (active_sites.size() == 0)
            return;

        BOOST_LOG_TRIVIAL(debug) << "Energy decreased, so assigning new labeling";
        for (VertexDescriptor vertex_desc : active_sites)
        {
            if (energy(vertex_desc).color == boost::black_color)
                continue;

            auto data_cost = 0;
            if (m_data_cost_fn) {
                auto vertex_idx = whichVertexIndex(vertex_desc);

                auto args = make_tuple(vertex_idx, alpha_label);
                data_cost = safeInvokeCostFn(m_data_cost_fn, args, 0);
            }

            m_sites_store.assignLabel(vertex_desc, alpha_label, data_cost);
        }
    }

    EnergyType computeEnergy()
    {
        EnergyType energy = 0;
        energy += m_runtime_statistics.pushEnergyToHistory("data", computeDataEnergy());
        energy += m_runtime_statistics.pushEnergyToHistory("smooth", computeSmoothEnergy());
        energy += m_runtime_statistics.pushEnergyToHistory("label", computeLabelEnergy());

        return energy;
    }

    EnergyType computeDataEnergy()
    {
        if (! m_data_cost_fn)
            return 0;

        EnergyType energy = 0;
        vector<VertexDescriptor> vertices;
        m_sites_store.queryAllVertices(vertices);

        for (auto vertex_desc : vertices)
            energy += m_sites_store.dataCost(vertex_desc);

        return energy;
    }

    EnergyType computeSmoothEnergy()
    {
        if (! m_smooth_cost_fn)
            return 0;

        return computeSmoothingTypeEnergy(m_smooth_cost_fn);
    }

    EnergyType computeLabelEnergy()
    {
        if (! m_label_cost_fn)
            return 0;

        return computeSmoothingTypeEnergy(m_label_cost_fn);
    }

    template<typename FnType>
    EnergyType computeSmoothingTypeEnergy(const FnType cost_fn)
    {
        using namespace std;

        EnergyType energy = 0;
        vector<VertexDescriptor> vertices;
        m_sites_store.queryAllVertices(vertices);

        for (auto vert_desc : vertices)
        {
            auto neighbouring_vertices = m_energy_graph.neighboursOf(vert_desc);
            for (auto nb_vert_desc : neighbouring_vertices)
            {
                auto cur_label = m_sites_store.whichLabel(vert_desc);
                auto nb_label = m_sites_store.whichLabel(nb_vert_desc);

                auto vert_idx = whichVertexIndex(vert_desc);
                auto nb_idx = whichVertexIndex(nb_vert_desc);

                auto args = make_tuple(vert_idx, nb_idx, cur_label, nb_label);
                energy += safeInvokeCostFn(cost_fn, args, 0);
            }
        }

        return energy;
    }

    inline int whichVertexIndex(const VertexDescriptor vert_desc)
    {
        auto node = m_energy_graph(vert_desc);
        auto idx = node.index;

        return idx - 2;
    }

    template<typename FnType, typename ArgType, typename CtxType>
    inline EnergyType safeInvokeCostFn(FnType fn, ArgType args, CtxType ctx)
    {
        EnergyType cost = fn(args, ctx);
        BOOST_ASSERT(cost >= 0);

        return cost;
    }

    void dumpEnergyGraph(int iter, int label_iter, int alpha_label, int energy, bool display = false)
    {
        if (! m_record_energy_graph_dumps)
            return;

        stringstream gvName;
        gvName << setw(3) << setfill('0') << iter << "_"
               << setw(3) << setfill('0') << label_iter << "_"
               << "label_" << setw(5) << setfill('0') << alpha_label << "_"
               << "energy_" << setw(5) << setfill('0') << energy << ".gv";
        m_energy_graph.dumpAsGraphviz(gvName.str());

        if (display)
        {
            string cmd = "dot -Tsvg " + gvName.str() + " | display";    
            assert(system(cmd.c_str()) == 0);
        }
    }

    void recordEnergyHistory(int iter, int label_iter, int alpha_label, int energy, bool display = false)
    {
        if (! m_record_energy_history)
            return;

        computeEnergy();
    }
};

#endif // BINOPT_H
