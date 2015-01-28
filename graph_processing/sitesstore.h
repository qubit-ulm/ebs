#ifndef SITESSTORE_H
#define SITESSTORE_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/density.hpp>

using namespace std;
namespace ba = boost::accumulators;

template<typename EnergyType = long long,
         typename VertexDescriptor = int,
         typename LabelType = int>
class SitesStore 
{
    struct Site 
	{
        VertexDescriptor m_vertex;
        LabelType m_label;
        bool m_is_active;
        EnergyType m_data_cost;
        EnergyType m_label_cost;

        Site(VertexDescriptor vertex_desc, LabelType label, bool is_active)
            : m_vertex(vertex_desc)
            , m_label(label)
            , m_is_active(is_active)
            , m_data_cost(0)
            , m_label_cost(0)
        { }

        Site(VertexDescriptor vertex_desc)
            : m_vertex(vertex_desc)
            , m_label(0)
            , m_is_active(false)
            , m_data_cost(0)
            , m_label_cost(0)
        { }
    };

    typedef boost::multi_index::multi_index_container<
      Site,
      boost::multi_index::indexed_by<
        boost::multi_index::hashed_unique<
          boost::multi_index::member<
            Site, VertexDescriptor, &Site::m_vertex
          >
        >,
        boost::multi_index::hashed_non_unique<
          boost::multi_index::member<
            Site, LabelType, &Site::m_label
          >
        >,
        boost::multi_index::hashed_non_unique<
            boost::multi_index::member<
                Site, bool, &Site::m_is_active
            >
        >
      >
    > SiteSet;
	
public:

    SitesStore()
		: m_vertex_index(m_sites.template get<0>())
		, m_label_index(m_sites.template get<1>())
		, m_active_index(m_sites.template get<2>())
    { }

    SitesStore& addVertices(vector<VertexDescriptor>& vertices)
    {
        for (auto vertex_desc : vertices)
            addVertex(vertex_desc);

        return *this;
    }

    SitesStore& addVertex(VertexDescriptor vertex_desc)
    {
        m_sites.insert(Site(vertex_desc));
        return *this;
    }

    bool hasVertex(VertexDescriptor vertex_desc)
    {
        auto it = m_vertex_index.find(vertex_desc);
        return it != m_vertex_index.end();
    }

    SitesStore& assignLabel(VertexDescriptor vertex_desc, LabelType label)
    {
        auto it = m_vertex_index.find(vertex_desc);
        m_vertex_index.modify(it, [label](Site& s) {
           s.m_label = label;
        });

        return *this;
    }

    SitesStore& assignDataCost(VertexDescriptor vertex_desc, EnergyType data_cost)
    {
        auto it = m_vertex_index.find(vertex_desc);
        m_vertex_index.modify(it, [data_cost](Site& s) {
           s.m_data_cost = data_cost;
        });

        return *this;
    }

    SitesStore& assignLabelCost(VertexDescriptor vertex_desc, EnergyType label_cost)
    {
        auto it = m_vertex_index.find(vertex_desc);
        m_vertex_index.modify(it, [label_cost](Site& s) {
           s.m_label_cost = label_cost;
        });

        return *this;
    }

    SitesStore& assignLabel(VertexDescriptor vertex_desc, LabelType label, long data_cost)
    {
        assignLabel(vertex_desc, label);
        assignDataCost(vertex_desc, data_cost);

        return *this;
    }

    SitesStore& setActive(VertexDescriptor vertex_desc, bool is_active)
    {
        auto it = m_vertex_index.find(vertex_desc);
        m_vertex_index.modify(it, [is_active](Site& s) {
           s.m_is_active = is_active;
        });

        return *this;
    }

    SitesStore& setActiveForLabel(LabelType label, bool is_active)
    {
        auto it = m_label_index.find(label);
        for (; it != m_label_index.end(); it++) {
            setActive((*it).m_vertex, is_active);
        }

        return *this;
    }

    LabelType whichLabel(VertexDescriptor vertex_desc) {
        auto it = m_vertex_index.find(vertex_desc);
        return (*it).m_label;
    }

    int dataCost(VertexDescriptor vertex_desc) {
        auto it = m_vertex_index.find(vertex_desc);
        return (*it).m_data_cost;
    }

    int labelCost(VertexDescriptor vertex_desc) {
        auto it = m_vertex_index.find(vertex_desc);
        return (*it).m_label_cost;
    }

    vector<VertexDescriptor> queryActiveVerticesForLabel(LabelType alpha_label)
    {
        vector<VertexDescriptor> active_vertices;

        auto site_iter = m_vertex_index.begin();
        for (; site_iter != m_vertex_index.end(); ++site_iter) {
            if ((*site_iter).m_label == alpha_label)
                continue;

            active_vertices.push_back((*site_iter).m_vertex);
        }

        return active_vertices;
    }

    void queryAllVertices(vector<VertexDescriptor>& vertices)
    {
        vertices.clear();

        auto site_iter = m_vertex_index.begin();
        for (; site_iter != m_vertex_index.end(); site_iter++) {
            vertices.push_back((*site_iter).m_vertex);
        }
    }

    SitesStore& markAllVerticesInactive() {
        auto site_iter = m_vertex_index.begin();
        for (; site_iter != m_vertex_index.end(); ++site_iter) {
            setActive((*site_iter).m_vertex, false);
        }

        return *this;
    }

    SitesStore& updateCountingStatistics()
    {
        updateLabelCounts();
        updateTransitionCount();

        return *this;
    }

    int labelCount(LabelType l)
    {
        if (! m_label_counts.count(l))
            return -1;

        return m_label_counts[l];
    }

    const std::map<LabelType, int>& getLabelCounts()
    {
        return m_label_counts;
    }

    const std::map<std::pair<LabelType, LabelType>, int>& getTransitionCounts()
    {
        return m_transition_counts;
    }



private:
    typedef typename SiteSet::template nth_index<0>::type SitesByVertex;
    typedef typename SiteSet::template nth_index<1>::type SitesByLabel;
    typedef typename SiteSet::template nth_index<2>::type SitesByActiveFlag;

    SiteSet m_sites;
    SitesByVertex& m_vertex_index;
    SitesByLabel& m_label_index;
    SitesByActiveFlag& m_active_index;

    map<LabelType, std::size_t> m_label_counts;
    map<pair<LabelType, LabelType>, int> m_transition_counts;

private:
    void updateLabelCounts()
    {
        using namespace std;

        m_label_counts.clear();

        for (auto vertex : m_vertex_index) {
            m_label_counts[vertex.m_label]++;
        }
    }

    void updateTransitionCount()
    {
        using namespace std;

        m_transition_counts.clear();

        shared_ptr<Site> prev_site;
        auto clique_count = 1;
        for (Site site : m_vertex_index) {
            if (! prev_site) {
                prev_site = make_shared<Site>(site);
                continue;
            }

            if (site.m_label == prev_site->m_label) {
                prev_site = make_shared<Site>(site);
                clique_count++;
                continue;
            }

            incrementTransitionCount(prev_site->m_label,
                                     site.m_label);

            clique_count = 1;
            prev_site = make_shared<Site>(site);
        }
    }

    inline void incrementTransitionCount(LabelType prev_label, LabelType cur_label)
    {
        auto key = make_pair(prev_label, cur_label);

        auto it = m_transition_counts.find(key);
        if (it == m_transition_counts.end()) {
            m_transition_counts[key] = 0;
        }

        m_transition_counts[key] += 1;
    }
};

#endif
