#ifndef ENERGY_H
#define ENERGY_H

#include <functional>

#include <boost/iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "helpers.h"
#include "bk_max_flow.h"

template<typename EnergyType = long long>
class EnergyGraph
{
public:
    typedef boost::mapS VertexContainerType;
    typedef boost::mapS EdgeContainerType;

    typedef boost::adjacency_list_traits<VertexContainerType, EdgeContainerType, boost::bidirectionalS> Traits;
    typedef typename Traits::vertex_descriptor VertexDescriptor;
    typedef typename Traits::edge_descriptor EdgeDescriptor;

    struct VertexProperties
    {
        VertexProperties()
            : index(-1)
            , color(boost::gray_color)
        { }

        VertexDescriptor desc;
        int index;
        std::string name;
        boost::default_color_type color;
    };

    struct EdgeProperties
    {
        EdgeProperties()
            : capacity(0)
            , residual_capacity(0)
        { }

        EdgeDescriptor reverse;
        EnergyType capacity;
        EnergyType residual_capacity;
    };

    struct GraphPropertiesType
    {
        EnergyType energy_const;
    };

    typedef boost::adjacency_list<VertexContainerType,
                                  EdgeContainerType,
                                  boost::bidirectionalS,
                                  VertexProperties,
                                  EdgeProperties,
                                  GraphPropertiesType
                                > Graph;

    typedef typename boost::property_map<Graph, VertexDescriptor VertexProperties::*>::type DescMapType;
    typedef typename boost::property_map<Graph, int VertexProperties::*>::type IndexMapType;
    typedef typename boost::property_map<Graph, std::string VertexProperties::*>::type NameMapType;
    typedef typename boost::property_map<Graph, boost::default_color_type VertexProperties::*>::type ColorMapType;
    typedef typename boost::property_map<Graph, EnergyType EdgeProperties::*>::type CapacityMapType;
    typedef typename boost::property_map<Graph, EdgeDescriptor EdgeProperties::*>::type ReverseMapType;
    typedef typename boost::property_map<Graph, EnergyType EdgeProperties::*>::type ResidualCapacityType;

    typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIter;
    typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIter;
    typedef typename boost::graph_traits<Graph>::adjacency_iterator AdjIter;

    enum EdgeDirection { IN, OUT };

public:
    EnergyGraph()
        : m_current_index(0)
        , m_energy_const(0)
        , m_flow(0)
        , m_check_submodularity(true)
    {
        initializePropertyMaps();
        initializeTerminalVertices();
    }

    ~EnergyGraph()
    { }

    inline std::pair<VertexIter, VertexIter> variableIterator()
    {
        return vertices(m_energy_graph);
    }

    inline std::vector<VertexDescriptor> neighboursOf(VertexDescriptor v)
    {
        std::vector<VertexDescriptor> neighbours;

        AdjIter it, it_end;
        boost::tie(it, it_end) = boost::adjacent_vertices(v, m_energy_graph);
        for(; it != it_end; it++) {
            if (*it == m_s_vertex || *it == m_t_vertex)
                continue;

            bool alreadyPresent = find(neighbours.begin(), neighbours.end(), *it) != neighbours.end();
            if (alreadyPresent)
                continue;

            neighbours.push_back(*it);
        }

        return neighbours;
    }

    inline VertexProperties& operator()(VertexDescriptor vert_desc)
    {
        return m_energy_graph[vert_desc];
    }

    inline VertexProperties& operator()(std::string name)
    {
        auto vert_decl = getVertexByName(name);
        return m_energy_graph[vert_decl];
    }

    long addConstant(EnergyType energy_to_add)
    {
        m_energy_const += energy_to_add;
        return m_energy_const;
    }

    VertexDescriptor sourceVertex() {
        return m_s_vertex;
    }

    VertexDescriptor targetVertex() {
        return m_t_vertex;
    }

    VertexDescriptor addVariable(std::string name)
    {
        auto vertex_desc = boost::add_vertex(m_energy_graph);
        initializeNewlyAddedVariable(vertex_desc, name);

        return vertex_desc;
    }

    VertexDescriptor addVariable()
    {
        auto vertex_desc = boost::add_vertex(m_energy_graph);
        initializeNewlyAddedVariable(vertex_desc, std::to_string(m_current_index - 2));

        return vertex_desc;
    }

    void addTerm1(VertexDescriptor vert_desc, EnergyType A, EnergyType B)
    {
        addTerminalCapacity(vert_desc, B, A);
    }

    void addTerm2(VertexDescriptor vert_u, VertexDescriptor vert_v,
                  EnergyType A, EnergyType B, EnergyType C, EnergyType D)
    {
        /* Proposed decomposition
        E = A B = A A  +  0  B-A
            C D   D D    C-D  0
        */

        /*
         A A
         D D
        */
        addTerminalCapacity(vert_u, 0, A);
        addTerminalCapacity(vert_v, D, 0);

        B = B - A;
        C = C - D;

        /* Remaining:
        0 B
        C 0
        */
        if (m_check_submodularity &&  B + C < 0)
        {
            std::stringstream msg;
            msg << "Suplied energy function is not regular and/or submodular. "
                << "B: " << B << ", C:" << C;

            throw new std::runtime_error(msg.str());
        }

        //assert(B + C >= 0); // regular, submodular function

        if (B < 0)
        {
            /*
            B B + -B 0 +  0   0
            0 0   -B 0   B+C  0
            */
            addTerminalCapacity(vert_u, 0,  B);
            addTerminalCapacity(vert_v, 0, -B);
            addEdge(vert_u, vert_v, 0, B+C);
        }
        else if (C < 0)
        {
            /*
            -C -C + C 0 +  0  B+C
             0  0   C 0    0   0
            */
            addTerminalCapacity(vert_u, 0, -C);
            addTerminalCapacity(vert_v, 0,  C);
            addEdge(vert_u, vert_v, B+C, 0);
        }
        else // B >= 0, C >= 0
        {
            addEdge(vert_u, vert_v, B, C);
        }
    }

    /*
    void addTerm3(VertexDescriptor u, VertexDescriptor v, VertexDescriptor w,
                  int E000, int E001,
                  int E010, int E011,
                  int E100, int E101,
                  int E110, int E111)
    {
        auto pi = (E000 + E011 + E101 + E110) - (E100 + E010 + E001 + E111);

        if (pi >= 0)
        {
            addConstant(E111 - (E011 + E101 + E110));
            addTerminalCapacity(u, E101, E001);
            addTerminalCapacity(v, E110, E100);
            addTerminalCapacity(w, E011, E010);

            auto delta = (E010 + E001) - (E000 + E011); // -pi(E[x=0])
            assert(delta >= 0); // check regularity
            addEdge(v, w, delta, 0);

            delta = (E100 + E001) - (E000 + E101); // -pi(E[y=0])
            assert(delta >= 0); // check regularity
            addEdge(w, u, delta, 0);

            delta = (E100 + E010) - (E000 + E110); // -pi(E[z=0])
            assert(delta >= 0); // check regularity
            addEdge(u, v, delta, 0);
        }
        else
        {
            addConstant(E000 - (E100 + E010 + E001));

            addTerminalCapacity(u, E110, E010);
            addTerminalCapacity(v, E011, E001);
            addTerminalCapacity(w, E101, E100);

            auto delta = (E110 + E101) - (E100 + E111); // -pi(E[x=1])
            assert(delta >= 0); // check regularity
            addEdge(w, v, delta, 0);

            delta = (E110 + E011) - (E010 + E111); // -pi(E[y=1])
            assert(delta >= 0); // check regularity
            addEdge(v, w, delta, 0);

            delta = (E101 + E011) - (E001 + E111); // -pi(E[z=1])
            assert(delta >= 0); // check regularity
            addEdge(v, u, delta, 0);
        }
    }
    */

    template <typename TIter = std::vector<EdgeDescriptor>::iterator>
    inline std::pair<TIter, TIter> termIterator(VertexDescriptor vert_desc)
    {
        std::vector<EdgeDescriptor> merged;

        auto out_iter = boost::out_edges(vert_desc, m_energy_graph);
        auto in_iter = boost::in_edges(vert_desc, m_energy_graph);

        std::merge(out_iter.first, out_iter.second,
                   in_iter.first, in_iter.second,
                   merged.begin());

        return std::make_pair(merged.begin(), merged.end());
    }

    inline EdgeProperties& operator()(VertexDescriptor vert_desc,
                                      EdgeDirection dir,
                                      int offset)
    {
        using namespace boost;

        if (dir == IN)
        {
            auto edge_iter = in_edges(vert_desc, m_energy_graph);
            return getNthEdge(edge_iter.first, edge_iter.second, offset);
        }
        else
        {
            auto edge_iter = out_edges(vert_desc, m_energy_graph);
            return getNthEdge(edge_iter.first, edge_iter.second, offset);
        }
    }

    template<typename TIterBegin, typename TIterEnd>
    inline EdgeProperties& getNthEdge(TIterBegin& iter_begin,
                                      TIterEnd& iter_end,
                                      int offset)
    {
        int i = 0;
        for (auto iter = iter_begin; iter != iter_end; iter++)
        {
            if (i == offset)
                return m_energy_graph[*iter];
            i++;
        }

        throw std::out_of_range("Unable to find offset");
    }

    inline EdgeProperties& operator()(VertexDescriptor src_vert_desc,
                                      EdgeDirection dir,
                                      VertexDescriptor tgt_vert_desc)
    {
        using namespace boost;

        if (dir == IN)
        {
            auto edge_iter = in_edges(src_vert_desc, m_energy_graph);
            return getFirstEdgeConnectedTo(edge_iter.first, edge_iter.second, tgt_vert_desc, dir);
        }
        else
        {
            auto edge_iter = out_edges(src_vert_desc, m_energy_graph);
            return getFirstEdgeConnectedTo(edge_iter.first, edge_iter.second, tgt_vert_desc, dir);
        }

    }

    template<typename TIterBegin, typename TIterEnd>
    inline EdgeProperties& getFirstEdgeConnectedTo(TIterBegin& iter_begin,
                                                   TIterEnd& iter_end,
                                                   VertexDescriptor tgt_vert_desc,
                                                   EdgeDirection dir)
    {
        for (auto iter = iter_begin; iter != iter_end; iter++)
        {
            auto cur_tgt_desc = dir == OUT
                                    ? boost::target(*iter, m_energy_graph)
                                    : boost::source(*iter, m_energy_graph);
            if (cur_tgt_desc == tgt_vert_desc)
                return m_energy_graph[*iter];
        }

        throw std::out_of_range("Unable to find an edge that connects to target vertex");
    }

    inline EdgeProperties& operator()(std::string vert_name,
                                          EdgeDirection dir,
                                          int offset)
    {
        auto vert_decl = getVertexByName(vert_name);
        return (*this)(vert_decl, dir, offset);
    }

    inline EdgeProperties& operator()(std::string src_vert_name,
                                          EdgeDirection dir,
                                          std::string tgt_vert_name)
    {
        auto src_vert_decl = getVertexByName(src_vert_name);
        auto tgt_vert_decl = getVertexByName(tgt_vert_name);

        return (*this)(src_vert_decl, dir, tgt_vert_decl);
    }

    EnergyType minimize()
    {
        using namespace boost;

        /*
        typedef std::map<VertexDescriptor, std::size_t> IndexMap;
        IndexMap map_index;

        boost::associative_property_map<IndexMap> propmap_index(map_index);
        //indexing the vertices
        int i=0;

        VertexIter vi, vi_end;
        tie(vi, vi_end) = vertices(m_energy_graph);

        for (; vi != vi_end; vi++)
        {
            boost::put(propmap_index, *vi, i++);
        }
        */
        auto max_flow_bk = bk_max_flow(m_energy_graph,
                                       m_capacity_prop,
                                       m_residual_capacity_prop,
                                       m_reverse_prop,
                                       m_color_prop,
                                       //get(boost::vertex_index, m_energy_graph),
                                       m_index_prop,
                                       m_s_vertex,
                                       m_t_vertex);

        return max_flow_bk + m_energy_const;
    }

    void recycle()
    {
        using namespace boost;

        VertexIter vi, vi_end;
        boost::tie(vi, vi_end) = vertices(m_energy_graph);

        EdgeIter ei, ei_end;
        boost::tie(ei, ei_end) = edges(m_energy_graph);

        for (; vi != vi_end; vi++)
            recycleVertex(*vi);

        for (; ei != ei_end; ei++)
            recycleEdge(*ei);
    }

    void dumpAsGraphviz(const std::string file_name)
    {
        std::ofstream ostream(file_name);
        dumpAsGraphviz(ostream);
        ostream.close();
    }

    void dumpAsGraphviz()
    {
        dumpAsGraphviz(std::cout);
    }

    void dumpAsGraphviz(std::ostream& ostream)
    {
        using namespace boost;

        default_writer dflt_writer;

        auto vertex_writer = [this] (std::ostream& os, VertexDescriptor u) {
            auto color = get(m_color_prop, u);

            os << " ["
               << "label=\"" << get(m_name_prop, u) << "\"";
               //<< "[" << get(m_desc_prop, u) << "]\"";
            if (color == black_color)
                os << ", color=black, fontcolor=white, style=filled";
            os << "]";
        };

        auto edge_writer = [this] (std::ostream& os, EdgeDescriptor e) {
            auto cap = get(m_capacity_prop, e);
            auto res_cap = get(m_residual_capacity_prop, e);
            auto flow = cap - res_cap;

            os << " [ "
               << "label=\"c:" << flow << "/" << cap << "\"";

            if (res_cap == 0 && cap == 0)
                os << ", style=invis";

            os << "]";
        };

        write_graphviz(ostream,
                       m_energy_graph,
                       vertex_writer,
                       edge_writer,
                       dflt_writer,
                       m_index_prop);
    }

private:
    Graph m_energy_graph;
    VertexDescriptor m_s_vertex;
    VertexDescriptor m_t_vertex;

    DescMapType m_desc_prop;
    int m_current_index;
    IndexMapType m_index_prop;
    NameMapType m_name_prop;
    ColorMapType m_color_prop;
    CapacityMapType m_capacity_prop;
    ReverseMapType m_reverse_prop;
    ResidualCapacityType m_residual_capacity_prop;

    EnergyType m_energy_const;
    EnergyType m_flow;

    bool m_check_submodularity;

private:
    void initializeTerminalVertices()
    {
        m_s_vertex = boost::add_vertex(m_energy_graph);
        put(m_desc_prop, m_s_vertex, m_s_vertex);
        put(m_index_prop, m_s_vertex, m_current_index++);
        put(m_name_prop, m_s_vertex, "s");

        m_t_vertex = boost::add_vertex(m_energy_graph);
        put(m_desc_prop, m_t_vertex, m_t_vertex);
        put(m_index_prop, m_t_vertex, m_current_index++);
        put(m_name_prop, m_t_vertex, "t");
    }

    void initializePropertyMaps()
    {
        using namespace boost;

        m_desc_prop = get(&VertexProperties::desc, m_energy_graph);
        m_index_prop = get(&VertexProperties::index, m_energy_graph);
        m_name_prop = get(&VertexProperties::name, m_energy_graph);
        m_color_prop = get(&VertexProperties::color, m_energy_graph);
        m_capacity_prop = get(&EdgeProperties::capacity, m_energy_graph);
        m_residual_capacity_prop = get(&EdgeProperties::residual_capacity, m_energy_graph);
        m_reverse_prop = get(&EdgeProperties::reverse, m_energy_graph);
    }

    void initializeNewlyAddedVariable(VertexDescriptor& vert_desc,
                                      std::string name,
                                      bool connect_to_terminals = true)
    {
        put(m_desc_prop, vert_desc, vert_desc);
        put(m_index_prop, vert_desc, m_current_index++);
        put(m_name_prop, vert_desc, name);
        put(m_color_prop, vert_desc, boost::white_color);

        if (connect_to_terminals) {
            addEdge(vert_desc, m_s_vertex, 0, 0);
            addEdge(vert_desc, m_t_vertex, 0, 0);
        }
    }

    inline void recycleVertex(const VertexDescriptor vert_desc)
    {
        using namespace boost;
        put(m_color_prop, vert_desc, boost::white_color);
    }

    inline void recycleEdge(const EdgeDescriptor edge_desc)
    {
        using namespace boost;

        put(m_capacity_prop, edge_desc, 0.0);
        put(m_residual_capacity_prop, edge_desc, 0.0);
    }

    VertexDescriptor getVertexByName(const std::string& name)
    {
        auto predicate = [this, name] (VertexDescriptor vert_desc) {
            auto vertex = m_energy_graph[vert_desc];
            return name == vertex.name;
        };

        return getVertexByPredicate(predicate);
    }

    VertexDescriptor getVertexByIndex(const int index)
    {
        auto predicate = [this, index] (VertexDescriptor vert_desc) {
            auto vertex = m_energy_graph[vert_desc];
            return index == vertex.index;
        };

        return getVertexByPredicate(predicate);
    }

    VertexDescriptor getVertexByPredicate(const std::function<bool(VertexDescriptor)>& predicate)
    {
        VertexIter vert_begin, vert_end;
        tie(vert_begin, vert_end) = variableIterator();

        auto filter_iter = boost::make_filter_iterator(predicate, vert_begin, vert_end);
        return *filter_iter;
    }

    inline void addTerminalCapacity(VertexDescriptor vertex_desc,
                                    long source_cap,
                                    long target_cap)
    {
        using namespace boost;

        EdgeDescriptor s_edge, t_edge;
        bool success = false;

        if (source_cap < 0) {
            target_cap += abs(source_cap);
            source_cap = 0;
        }

        if (target_cap < 0) {
            source_cap += abs(target_cap);
            target_cap = 0;
        }

        boost::tie(s_edge, success) = boost::edge(m_s_vertex, vertex_desc, m_energy_graph);
        if (! success && source_cap != 0) {
            boost::tie(s_edge, success) = add_edge(m_s_vertex, vertex_desc, m_energy_graph);
            setEdgeCapacity(s_edge, 0);
        }

        boost::tie(t_edge, success) = boost::edge(vertex_desc, m_t_vertex, m_energy_graph);
        if (! success && target_cap != 0) {
            boost::tie(t_edge, success) = add_edge(vertex_desc, m_t_vertex, m_energy_graph);
            setEdgeCapacity(t_edge, 0);
        }

        if (source_cap != 0)
            setEdgeCapacity(s_edge, source_cap);

        if (target_cap != 0)
            setEdgeCapacity(t_edge, target_cap);
    }

    void addEdge(const VertexDescriptor& vert_u, const VertexDescriptor& vert_v,
                 long cap, long rev_cap)
    {
        using namespace boost;

        EdgeDescriptor edge_uv, edge_vu;
        bool success;

        if (cap < 0) {
            rev_cap += abs(cap);
            cap = 0;
        }

        if (rev_cap < 0) {
            cap += abs(rev_cap);
            rev_cap = 0;
        }

        boost::tie(edge_uv, success) = edge(vert_u, vert_v, m_energy_graph);
        if (! success) {
            boost::tie(edge_uv, success) = add_edge(vert_u, vert_v, m_energy_graph);
            setEdgeCapacity(edge_uv, 0);
        }

        boost::tie(edge_vu, success) = edge(vert_v, vert_u, m_energy_graph);
        if (! success) {
            boost::tie(edge_vu, success) = add_edge(vert_v, vert_u, m_energy_graph);
            setEdgeCapacity(edge_vu, 0);
        }

        setEdgeCapacity(edge_uv, cap);
        setReverseEdge(edge_uv, edge_vu);

        setEdgeCapacity(edge_vu, rev_cap);
        setReverseEdge(edge_vu, edge_uv);
    }

    inline void setEdgeCapacity(const EdgeDescriptor& edge_desc, long value)
    {
        auto delta = get(m_capacity_prop, edge_desc);
        value += delta;

        put(m_capacity_prop, edge_desc, value);
    }

    inline void setReverseEdge(const EdgeDescriptor& edge_desc, const EdgeDescriptor &rev_edge)
    {
        put(m_reverse_prop, edge_desc, rev_edge);
    }
};

#endif // ENERGY_H
