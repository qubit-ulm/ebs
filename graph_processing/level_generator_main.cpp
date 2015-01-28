#include <algorithm>
#include <iostream>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "cmd_helpers.h"

namespace bpo = boost::program_options;

namespace 
{
    bpo::options_description initializeOptionsDescription()
    {
        using namespace std;

        bpo::options_description desc("Options");
        desc.add_options()
            ("help,h", "Print help message")
            ("input", bpo::value<string>(),
                    "Filename of a matrix market vector file "
                    "containing the denoised input data set")
            ("output", bpo::value<string>()->default_value("-"),
                    "Filename of the matrix market vector file "
                    "the level data should be written to")
            ("level-distance", bpo::value<double>(),
                    "Distance between each level between "
                    "min/max value of input vector")
            ("level-number", bpo::value<size_t>(),
                    "Number of linearly spaced leves between "
                    "min/max value of input vector")
            ("debug,d", "Turn on debug output if flag is set");

        return desc;
    }

    bool tryParseProgramOptions(const bpo::options_description &desc,
                                const int argc,
                                const char *argv[],
                                bpo::variables_map &vm)
    {
        using namespace std;
        using namespace boost::program_options;
        
        bpo::positional_options_description ppos;
        ppos.add("input", 1);
        ppos.add("output", 2);

        try {
            bpo::store(bpo::command_line_parser(argc, argv)
                        .options(desc)
                        .positional(ppos)
                        .run(), vm);
            bpo::notify(vm);


            if (cmd::isHelpRequest(desc, vm)) {
                return cmd::SUCCESS;
            }
            
            bool is_valid = true;

            if (! vm.count("input")) {
                cout << "ERROR: 'input' argument is required" << endl;
                is_valid = false;
            }

            if (! vm.count("level-distance") && ! vm.count("level-number")) {
                cout << "ERROR: either 'level-distance' or 'level-number' argument is required" << endl;
                is_valid = false;
            }

            if (! vm.count("output")) {
                cout << "ERROR: 'output' argument is required" << endl;
                is_valid = false;
            }

            return is_valid;
        }
        catch (error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            return false;
        }

        return false;
    }

    template <typename VectorType>
    std::pair<double, double> getVectorMinMax(const VectorType &v)
    {
        using namespace std;

        auto min = min_element(v.begin(), v.end());
        auto max = max_element(v.begin(), v.end());
        return make_pair(*min, *max);
    }

    template <typename VectorType>
    void populateLevels(VectorType &levels,
                        double min, 
                        double distance,
                        size_t n)
    {
        levels.resize(n);

        for (size_t i = 0; i < n; i++)
            levels(i) = min + distance * i;
    }

    template <typename VectorType>
    void populateLevelsByDistance(const bpo::variables_map &vm,
                                  const VectorType &input,
                                  VectorType &levels)
    {
        using namespace std;

        double distance = vm["level-distance"].as<double>();

        double min, max;
        tie(min, max) = getVectorMinMax(input);

        size_t n = ceil((max - min) / distance);

        populateLevels(levels, min, distance, n);
    }

    template <typename VectorType>
    void populateLevelsByNumber(const bpo::variables_map &vm,
                                const VectorType &input,
                                VectorType &levels)
    {
        using namespace std;

        auto n = vm["level-number"].as<size_t>();

        double min, max;
        tie(min, max) = getVectorMinMax(input);

        double distance = (max - min) / (double)n;

        populateLevels(levels, min, distance, n);
    }
                                

    int runProgram(const bpo::variables_map &vm)
    {
        using namespace std;
        typedef boost::numeric::ublas::vector<double> VectorType;

        VectorType input;

        if (! cmd::loadInputVector(vm, input))
            return cmd::ERROR_UNHANDLED_EXCEPTION;

        VectorType levels;
        if (vm.count("level-distance"))
            populateLevelsByDistance(vm, input, levels);
        else
            populateLevelsByNumber(vm, input, levels);

        if (! cmd::saveOutputVector(vm, levels))
            return cmd::ERROR_UNHANDLED_EXCEPTION;

        return cmd::SUCCESS;
    }
}

int main (const int argc, const char *argv[])
{
    using namespace std;

    auto desc = initializeOptionsDescription();

    try 
    {
        bpo::variables_map vm;
        if (! tryParseProgramOptions(desc, argc, argv, vm)) {
            return cmd::ERROR_IN_COMMAND_LINE;
        }

        return runProgram(vm);
    }
    catch (exception &e)
    {
        cerr << "Unhandled Exception reached the top of main: "
             << e.what() << ", application will now exit" << endl;
    }
}
