#ifndef RUNTIME_STATISTICS_H
#define RUNTIME_STATISTICS_H

#include <map>
#include <vector>

template<typename EnergyLabel = std::string,
         typename EnergyType = long long>
class RuntimeStatistics
{
public:

    typedef std::map<EnergyLabel, std::vector<EnergyType>> EnergyHistoryType;

    RuntimeStatistics()
    { }

    EnergyType pushEnergyToHistory(EnergyLabel label, EnergyType energy)
    {
        ensureLabelVector(label);
        m_energy_history[label].push_back(energy);

        return energy;
    }

    std::vector<EnergyType>& findEnergyHistory(EnergyLabel label)
    {
        if (! energyHistoryContainsKey(label))
            throw std::out_of_range("Energy label does not exist");

        return m_energy_history.find(label)->second;
    }

private:
    EnergyHistoryType m_energy_history;

private:

    inline bool energyHistoryContainsKey(EnergyLabel label)
    {
        return m_energy_history.find(label) != m_energy_history.end();
    }

    void ensureLabelVector(EnergyLabel label)
    {
        if (energyHistoryContainsKey(label))
            return;

        m_energy_history[label] = std::vector<EnergyType>();
    }
};


#endif // RUNTIME_STATISTICS_H
