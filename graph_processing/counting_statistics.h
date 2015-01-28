#ifndef COUNTING_STATISTICS_H
#define COUNTING_STATISTICS_H

#include <cfloat>
#include <functional>
#include <random>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/density.hpp>

#include <boost/log/trivial.hpp>

#define BOOST_RESULT_OF_USE_DECLTYPE
#include <boost/iterator.hpp>
#include <boost/range.hpp>

namespace detail {

    template <typename TRet>
    TRet meanHelper(long double sum, std::size_t count)
    {
        if (count == 0)
            return TRet(0);

        const long double countf = count;
        return static_cast<TRet>(sum / countf);
    }

    template<typename TValue>
    struct Bin
    {
    public:
        Bin()
            : sum(TValue())
            , count(0)
            , weight(0.0)
        {}

        void clear()
        {
            sum = TValue();
            count = 0;
        }

        void add(const TValue& s, std::size_t c)
        {
            sum += s;
            count += c;
        }

        Bin& operator+=(const Bin& o)
        {
            add(o.sum, o.count);
            return *this;
        }

        Bin& operator-=(const Bin& o)
        {
            sum -= o.sum;
            count -= o.count;
            return *this;
        }

        TValue mean() const
        {
            return meanHelper<TValue>(sum, count);
        }

        TValue sum;
        std::size_t count;
        double weight;
    };

    template<typename TValue, typename TBin>
    class HistogramBinContainer
    {
    public:
        HistogramBinContainer(TValue bin_size,
                              TValue min_value,
                              TValue max_value,
                              const TBin& default_bin)
            : m_bin_size(bin_size)
            , m_min(min_value)
            , m_max(max_value)
            , m_event_count(0)
        {
            unsigned int number_bins = (max_value - min_value) / bin_size;

            // Round up if the bucket size does not fit evenly
            if (number_bins * bin_size < max_value - min_value)
                ++number_bins;

            // Add 2 for the extra 'below min' and 'above max' buckets
            number_bins += 2;
            m_bins.assign(number_bins, default_bin);
        }

        TValue getBinSize() const
        {
            return m_bin_size;
        }

        TValue getMin() const
        {
            return m_min;
        }

        TValue getMax() const
        {
            return m_max;
        }

        std::size_t getNumBins() const
        {
            return m_bins.size();
        }

        std::size_t getBinIdx(TValue value)
        {
            if (value < m_min)
                return 0;

            if (value >= m_max)
                return m_bins.size() - 1;

            return ((value - m_min) / m_bin_size) + 1;
        }

        TBin& getByValue(TValue value)
        {
            return m_bins[getBinIdx(value)];
        }

        TBin& getByIndex(std::size_t  idx)
        {
            auto& bin = m_bins[idx];
            updateBinWeight(bin);

            return bin;
        }

        TBin& operator[] (const std::size_t
                          idx)
        {
            return getByIndex(idx);
        }

        TValue getBinLowerBound(std::size_t idx) const
        {
            if (idx == 0)
                return std::numeric_limits<TValue>::min();

            if (idx == m_bins.size() - 1)
                return m_max;

            return m_min + ((idx - 1) * m_bin_size);
        }

        TValue getBinUpperBound(std::size_t idx) const
        {
            if (idx == m_bins.size() - 1)
                return std::numeric_limits<TValue>::max();

            return m_min + (idx * m_bin_size);
        }

        TValue getBinMeanValue(std::size_t idx) const
        {
            auto lower_bound = getBinLowerBound(idx);
            auto upper_bound = getBinUpperBound(idx);

            auto mean = lower_bound + (upper_bound - lower_bound) / 2;
            return min(getMax(), max(getMin(), mean));
        }

        std::size_t incEventCount(int n = 1) {
            m_event_count += n;
            return getNumEvents();
        }

        std::size_t decEventCount(int n = 1) {
            m_event_count -= n;
            if (m_event_count < 0)
                m_event_count = 0;

            return getNumEvents();
        }

        inline std::size_t getNumEvents() const {
            return m_event_count;
        }

        typedef std::function<TBin& (TBin&)> WeightTransformer;
        typedef boost::transform_iterator<WeightTransformer, typename std::vector<TBin>::iterator> iterator;

        HistogramBinContainer<TValue, TBin>::iterator begin()
        {
            WeightTransformer transformer = std::bind(&HistogramBinContainer::updateBinWeight, this, std::placeholders::_1);
            return boost::make_transform_iterator(m_bins.begin(), transformer);
        }

        HistogramBinContainer<TValue, TBin>::iterator end()
        {
            WeightTransformer transformer = std::bind(&HistogramBinContainer::updateBinWeight, this, std::placeholders::_1);
            return boost::make_transform_iterator(m_bins.end(), transformer);
        }

    private:
        TValue m_bin_size;
        TValue m_min;
        TValue m_max;
        std::vector<TBin> m_bins;

        std::size_t m_event_count;

    private:
        TBin& updateBinWeight(TBin& bin) const
        {
            bin.weight = (double)(bin.count) / getNumEvents();
            return bin;
        }
    };
}

template<typename TValue = double>
class Histogram
{
    typedef detail::Bin<TValue> Bin;

public:
    Histogram(TValue bin_size,
              TValue min,
              TValue max)
        : m_bins(bin_size, min, max, detail::Bin<TValue>())
    { }

    void addValue(TValue value)
    {
        Bin& bin = m_bins.getByValue(value);
        bin.sum += value;
        bin.count += 1;

        m_bins.incEventCount();
    }

    void addRepeatedValue(TValue value, std::size_t n_samples)
    {
        Bin& bin = m_bins.getByValue(value);
        bin.sum += value * n_samples;
        bin.count += n_samples;

        m_bins.incEventCount(n_samples);
    }

    void removeValue(TValue value)
    {
        Bin& bin = m_bins.getByValue(value);

        if (bin.count > 0) {
            bin.sum -= value;
            bin.count -= 1;
        }
        else {
            bin.sum = TValue();
            bin.count = 0;
        }

        m_bins.decEventCount();
    }

    void removeRepeatedValue(TValue value, std::size_t n_samples)
    {
        Bin& bin = m_bins.getByValue(value);

        if (bin.count >= n_samples) {
            bin.sum -= value * n_samples;
            bin.count -= n_samples;
        }
        else {
            bin.sum = TValue();
            bin.count = 0;
        }

        m_bins.decEventCount(n_samples);
    }

    void clear()
    {
        for (auto i = 0; i < m_bins.getNumBins(); i++) {
            m_bins.getByIndex(i)
                  .clear();
        }

        m_bins.decEventCount(std::numeric_limits<std::size_t>::max());
    }

    bool doHistogramDimensionsMatch(const Histogram &hist)
    {
        return getBinSize() == hist.getBinSize() &&
               getMin() == hist.getMin() &&
               getMax() == hist.getMax() &&
               getNumBins() == hist.getNumBins();
    }

    void subtract(const Histogram &hist)
    {
        if (! doHistogramDimensionsMatch(hist))
            throw std::invalid_argument("Cannot subtract input histogram.");

        for (auto i = 0; i < getNumBins(); i++)
            getBinByIndex(i) -= hist.getBinByIndex(i);
    }

    void merge(const Histogram &hist)
    {
        if (! doHistogramDimensionsMatch(hist))
            throw std::invalid_argument("Cannot merge from input histogram.");

        for (auto i = 0; i < getNumBins(); i++)
            getBinByIndex(i) += hist.getBinByIndex(i);
    }

    void copy(const Histogram &hist)
    {
        if (! doHistogramDimensionsMatch(hist))
            throw std::invalid_argument("Cannot copy from input histogram.");

        for (int i = 0; i < getNumBins(); i++)
            getBinByIndex(i) = hist.getBinByIndex(i);
    }

    TValue getBinSize() const
    {
        return m_bins.getBinSize();
    }

    TValue getMin() const
    {
        return m_bins.getMin();
    }

    TValue getMax() const
    {
        return m_bins.getMax();
    }

    TValue getMean()
    {
        TValue mean = 0.0;
        for (auto i = 0; i < m_bins.getNumBins(); i++)
            mean += m_bins[i].weight * m_bins.getBinMeanValue(i);

        return mean;
    }

    TValue getVariance()
    {
        TValue hist_mean = getMean();

        TValue var = 0.0;
        for (auto i = 0; i < m_bins.getNumBins(); i++) {
            var += m_bins[i].weight * pow(m_bins.getBinMeanValue(i) - hist_mean, 2.0);
        }

        return sqrt(var - 1.0 / 12.0 * pow(m_bins.getBinSize(), 2));
    }

    std::size_t getNumBins() const
    {
        return m_bins.getNumBins();
    }

    detail::Bin<TValue>& getBinByValue(TValue value)
    {
        auto idx = m_bins.getBinIdx(value);
        return getBinByIndex(idx);
    }

    detail::Bin<TValue>& getBinByIndex(std::size_t idx)
    {
        detail::Bin<TValue>& bin = m_bins.getByIndex(idx);
        return bin;
    }

    detail::Bin<TValue>& operator[] (const int idx)
    {
        return getBinByIndex(idx);
    }

    TValue getBinLowerBound(std::size_t idx) const
    {
        return m_bins.getBinLowerBound(idx);
    }

    TValue getBinUpperBound(std::size_t idx) const
    {
        return m_bins.getBinUpperBound(idx);
    }

    std::size_t getNumEvents() const
    {
        return m_bins.getNumEvents();
    }

    typename detail::HistogramBinContainer<TValue, Bin>::iterator begin()
    {
        return m_bins.begin();
    }

    typename detail::HistogramBinContainer<TValue, Bin>::iterator end()
    {
        return m_bins.end();
    }

    enum CompareType { CHI_SQUARE, CHI_SQUARE_ALT, CORELL, BHATTACHARYYA };

    double compareTo(Histogram<TValue>& h2, CompareType compare_type = CHI_SQUARE)
    {
        switch (compare_type)
        {
            case CHI_SQUARE_ALT:
                return compareHistChiSquareAlt(*this, h2);
            case CORELL:
                return compareHistCorrel(*this, h2);
            case BHATTACHARYYA:
                return compareHistBhattacharyya(*this, h2);
            default:
                return compareHistChiSquare(*this, h2);
        }
    }

    void debugToStream(std::ostream& stream)
    {
        using namespace std;

        for (auto i = 0; i < getNumBins(); i++) {
            auto low = getBinLowerBound(i);
            auto up = getBinUpperBound(i);
            auto weight = max(0, (int)floor(getBinByIndex(i).weight * 100));

            stream << "[" << setw(8) << setprecision(2) << low
                   << "|" << setw(8) << setprecision(2) << up
                   << "] "
                   << string(weight, '*')
                   << endl;
        }
    }

private:
    detail::HistogramBinContainer<TValue, Bin> m_bins;

private:
    double compareHistChiSquare(Histogram<TValue>& h1, Histogram<TValue>& h2)
    {
        double result = 0.0;
        assert(h1.getNumBins() == h2.getNumBins());

        for (auto i = 0; i < h1.getNumBins(); i++)
        {
            double w1 = h1[i].weight;
            double w2 = h2[i].weight;

            double a = w1 - w2;
            double b = w1;

            if (fabs(b) > DBL_EPSILON)
                result += a * a / b;
        }

        return result;
    }

    double compareHistChiSquareAlt(Histogram<TValue>& h1, Histogram<TValue>& h2)
    {
        double result = 0.0;
        assert(h1.getNumBins() == h2.getNumBins());

        for (auto i = 0; i < h1.getNumBins(); i++)
        {
            double w1 = h1[i].weight;
            double w2 = h2[i].weight;

            double a = w1 - w2;
            double b = w1 + w2;

            if (fabs(b) > DBL_EPSILON)
                result += a * a / b;
        }

        return (result * 2);
    }

    double compareHistCorrel(Histogram<TValue>& h1, Histogram<TValue>& h2)
    {
        double result = 0.0;
        assert(h1.getNumBins() == h2.getNumBins());

        double s1 = 0, s2 = 0, s11 = 0, s12 = 0, s22 = 0;

        for (auto i = 0; i < h1.getNumBins(); i++)
        {
            double w1 = h1[i].weight;
            double w2 = h2[i].weight;

            double a = w1;
            double b = w2;

            s12 += a*b;
            s1 += a;
            s11 += a*a;
            s2 += b;
            s22 += b*b;
        }

        size_t total = h1.getNumBins();
        double scale = 1./total;
        double num = s12 - s1*s2*scale;
        double denom2 = (s11 - s1*s1*scale)*(s22 - s2*s2*scale);

        result = std::abs(denom2) > DBL_EPSILON
                ? num/std::sqrt(denom2)
                : 1.;

        return result;
    }

    double compareHistBhattacharyya(Histogram<TValue>& h1, Histogram<TValue>& h2)
    {
        double result = 0.0;
        assert(h1.getNumBins() == h2.getNumBins());

        double s1 = 0, s2 = 0;

        for (auto i = 0; i < h1.getNumBins(); i++)
        {
            double w1 = h1[i].weight;
            double w2 = h2[i].weight;

            double a = w1;
            double b = w2;

            result += std::sqrt(a*b);
            s1 += a;
            s2 += b;
        }

        s1 *= s2;
        s1 = fabs(s1) > DBL_EPSILON ? 1./std::sqrt(s1) : 1.;
        result = std::sqrt(std::max(1. - result*s1, 0.));

        return result;
    }
};

template<typename LabelDesc = int,
         typename LabelType = double>
class CountingStatistics
{
public:
    typedef std::map<std::pair<LabelDesc, LabelDesc>, int> TransitionHistoType;

    CountingStatistics(const double mean,
                       const double variance,
                       const std::vector<LabelType>& labels)
        : m_reference_mean(mean)
        , m_reference_variance(variance)
        , m_labels(labels)
        , m_bin_count(100)
        , m_reference_jump_histogram(createReferenceJumpHistogram())
    {
    }

    double penalizeTransitionConfiguration(const TransitionHistoType& transition_histogram)
    {
        using namespace std;

        double min_cost = 0.0;
        double max_cost = 100000.0;

        /*
        cout << setw(11) << "Trans.: "
             << setw(8) << setprecision(2) << transition_histogram.size()
             << endl;
        */

        auto jump_histogram = convertTransitionsToJumps(transition_histogram);
        //debugHistogram(jump_histogram, std::cout, false);

        /*
        if (jump_histogram.getNumEvents() < 15) {
            BOOST_LOG_TRIVIAL(debug) << "\tNum events below thres. ("
                                     <<  jump_histogram.getNumEvents()<< ")";
            return max_cost;
        }
        */

        auto cost = m_reference_jump_histogram.compareTo(jump_histogram);
        if (std::isnan(cost)) {
            BOOST_LOG_TRIVIAL(debug) << "\tLabel cost is NaN";
            return max_cost;
        }

        //BOOST_LOG_TRIVIAL(debug) << "\tLabel cost: "
        //       << setw(8) << setprecision(2) << cost;

        //return min(abs(cost), max_cost);
        return abs(cost);
    }

    void debugHistogram(Histogram<double>& histogram,
                        std::ostream& stream,
                        bool draw_ascii_graph = false)
    {
        using namespace std;

        if (draw_ascii_graph)
            histogram.debugToStream(stream);

         BOOST_LOG_TRIVIAL(debug)
               << "NEvt: "
               << setw(8) << setprecision(2) << histogram.getNumEvents()
               << "\tMean: "
               << setw(8) << setprecision(2) << histogram.getMean()
               << " ("
               << setw(8) << setprecision(2) << m_reference_jump_histogram.getMean()
               << ")"
               << "\tVar:  "
               << setw(8) << setprecision(2) << histogram.getVariance()
               << " ("
               << setw(8) << setprecision(2) << m_reference_jump_histogram.getVariance()
               << ")";
    }

    void debugReferenceHistogram(std::ostream& stream) {
       debugHistogram(m_reference_jump_histogram, stream, true);
    }

private:
    const std::vector<LabelType>& m_labels;
    const std::size_t m_bin_count;

    const double m_reference_mean;
    const double m_reference_variance;
    Histogram<double> m_reference_jump_histogram;

private:

    Histogram<double> createReferenceJumpHistogram()
    {
        using namespace std;

        auto jump_histogram = createHistogramConfiguredByLabel();

        random_device rd;
        mt19937 gen(rd());
        normal_distribution<double> dist(m_reference_mean, m_reference_variance);

        for (int i = 0; i < 10000; i++)
            jump_histogram.addValue(dist(gen));

        return jump_histogram;
    }

    Histogram<double> convertTransitionsToJumps(const TransitionHistoType& transition_histogram)
    {
        using namespace std;

        auto jump_histogram = createHistogramConfiguredByLabel();
        for(auto kv : transition_histogram)
        {
            auto trans = kv.first;
            auto site_count = kv.second;

            auto jump_height = m_labels[get<1>(trans)] - m_labels[get<0>(trans)];
            jump_histogram.addRepeatedValue(jump_height, site_count);
        }

        return jump_histogram;
    }

    std::pair<LabelType, LabelType> getLabelBounds()
    {
        auto minmax_pair = std::minmax_element(m_labels.begin(), m_labels.end());
        return make_pair(*(minmax_pair.first), *(minmax_pair.second));
    }

    Histogram<double> createHistogramConfiguredByLabel()
    {
        double bin_size, min_lbl, max_lbl;
        std::tie(min_lbl, max_lbl) = getLabelBounds();
        bin_size = (max_lbl - min_lbl) / m_bin_count;

        return Histogram<double>(bin_size, 0, max_lbl - min_lbl);
    }
};


#endif
