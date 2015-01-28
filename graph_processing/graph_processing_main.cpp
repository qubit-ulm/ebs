#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "../common/tuple_helper.h"
#include "cmd_helpers.h"
#include "binopt.h"

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
            ("levels", bpo::value<string>(),
                "Filename of a matrix market vector file "
                "containing the level set to culster the datapoints to")
            ("output", bpo::value<string>()->default_value("-"),
                    "Filename of the matrix market vector file "
                    "the denoised data should be written to")
            ("rho-d", bpo::value<double>()->default_value(100.0),
                    "Value of the regularization parameter for the data term")
            ("rho-s", bpo::value<double>()->default_value(10.0),
                    "Value of the regularization parameter for the smoothing term")
            ("rho-p", bpo::value<double>()->default_value(0.0),
                    "Value of the regularization parameter for the prior term")
            ("assignments", "Output assignments to levels and not the whole vector")
            ("maxiter", bpo::value<int>()->default_value(-1),
                    "Number of alpha expansion iterations, if set to -1 (default) "
                    "a backtracking level proposal strategy is used.")
            ("prior-distance", bpo::value<double>(),
                    "The distance of two adjacent steps the prior term should "
                    "NOT penalize")
            ("debug,d", "Turn on debug output if flag is set");
            ("debug-graphstructure", 
                    "Turn on debug output of graphs and capacities after "
                    "each graph-cut. The files are in graphviz dot notation.");

        return desc;
    }


    bool tryParseProgramOptions(const bpo::options_description& desc,
                                const int argc,
                                const char *argv[],
                                bpo::variables_map& vm)
    {
        using namespace std;
        using namespace boost::program_options;

        bpo::positional_options_description ppos;
        ppos.add("input", 1);
        ppos.add("levels", 2);
        ppos.add("output", 3);

        try {
            store(parse_command_line(argc, argv, desc), vm); // can throw

            if (cmd::isHelpRequest(desc, vm)) {
                return cmd::SUCCESS;
            }

            bool is_valid = true;
            if (! vm.count("input")) {
                cout << "ERROR: 'input' argument is required" << endl;
                is_valid = false;
            }

            if (! vm.count("levels")) {
                cout << "ERROR: 'levels' argument is required" << endl;
                is_valid = false;
            }

            if (! vm.count("output")) {
                cout << "ERROR: 'output' argument is required" << endl;
                is_valid = false;
            }

            if (vm["rho-p"].as<double>() != 0.0 && ! vm.count("prior-distance")) {
                cout << "ERROR: 'prior-distance' is required, if prior term is activated ";
                cout << "('rho-p' > 0.0)" << endl;
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
    }

    template <typename VectorType>
    void registerCostFunctions(BinaryOptimization<> &bin_opt,
                               VectorType &data,
                               VectorType &weights,
                               VectorType &labels,
                               VectorType &lambdas,
                               double prior_distance)
    {
        typedef std::tuple<int, int> DataTuple;
        typedef std::tuple<int, int, int, int> SmoothTuple;

        auto data_cost_fn = [&data, &weights, &labels, &lambdas](const DataTuple &t, int)
        {
            double value  = data(get<0>(t));
            double weight = weights(get<0>(t));
            double label_value = labels(get<1>(t));

            double cost  = lambdas[0]
                * (1.0 + weight)
                * abs(value - label_value);

            BOOST_ASSERT(cost >= 0);
            return (long long)(cost);
        };

        auto smooth_cost_fn = [&data, &weights, &labels, &lambdas](const SmoothTuple &t, int)
        {
            double weight_1  = weights(get<0>(t));
            double weight_2  = weights(get<1>(t));
            double label_1   = labels(get<2>(t));
            double label_2   = labels(get<3>(t));

            double cost = lambdas[1]
                * (1.0 + weight_1 + weight_2)
                * (label_1 != label_2 ? 1 : 0);

            BOOST_ASSERT(cost >= 0);
            return (long long)(cost);
        };

        auto label_cost_fn = [&data, &weights, &labels, &prior_distance, &lambdas](const SmoothTuple &t, int)
        {
            double weight_1     = weights(get<0>(t));
            double weight_2     = weights(get<1>(t));
            double label_1      = labels(get<2>(t));
            double label_2      = labels(get<3>(t));
            double label_value1 = labels(label_1);
            double label_value2 = labels(label_2);

            double epsilon = 0.05;
            double delta = abs(prior_distance - abs(label_value1 - label_value2));

            double cost = lambdas[2]
                * (delta > epsilon ? 1 : 0)*(label_1 != label_2 ? 1 : 0);

            BOOST_ASSERT(cost >= 0);
            return (long long)(cost);
        };

        bin_opt.setDataCost(data_cost_fn);
            bin_opt.setSmoothnessCost(smooth_cost_fn);

        if (fabs(lambdas[2]) > std::numeric_limits<double>::epsilon())
            bin_opt.setLabelCost(label_cost_fn);
    }

    bool areAssignmentsRequested(const bpo::variables_map &vm)
    {
        if (vm.count("assignments"))
            return true;

        return false;
    }

    template <typename MatType, typename VectorType>
    void collectAssignments(BinaryOptimization<> &bin_opt,
                            const VectorType &weights,
                            const VectorType &labels,
		    	            MatType &mat)
    {
        using namespace std;

        auto assignments = bin_opt.whichLabels();
        mat.resize(assignments.size(), 2);
        
        for (size_t row = 0; row < assignments.size(); row++)
        {
            size_t lbl_idx = assignments[row];

            mat(row, 0) = labels(lbl_idx);
            mat(row, 1) = weights(row);
        }
    }

    template <typename VectorType>
    bool saveAssignments(const bpo::variables_map &vm,
                         const VectorType &weights,
                         const VectorType &labels,
		                 BinaryOptimization<> &bin_opt)
    {
        using namespace std;

        boost::numeric::ublas::matrix<double> assignments;
        collectAssignments(bin_opt, weights, labels, assignments);

        return cmd::saveOutputMatrix(vm, assignments); 
    }

    template <typename VectorType>
    void postprocessAssignments(const bpo::variables_map &vm,
		    		            BinaryOptimization<> &bin_opt,
				                const VectorType &input,
				                const VectorType &weights,
				                const VectorType &labels,
				                VectorType &output)
    {
        using namespace std;

        output.resize(input.size());
        
        boost::numeric::ublas::matrix<double> assignments;
        collectAssignments(bin_opt, weights, labels, assignments);

        size_t i = 0;
        double total_weight = 0.0;
        for (size_t j = 0; j < assignments.size1(); j++)
        {
            double val = assignments(j, 0);
            double weight = assignments(j, 1);

            total_weight += weight;

            for (size_t w = 0; w < weight; w++){
                output(i) = val;
                i++;
            }
        }

        //BOOST_LOG_TRIVIAL(debug) << "Total weight: " << total_weight;
        //BOOST_LOG_TRIVIAL(debug) << "Output size:  " << output.size();
    }

    int runProgram(const bpo::options_description& desc,
                   const bpo::variables_map& vm)
    {
        using namespace std;
        typedef boost::numeric::ublas::vector<double> VectorType;
        
        VectorType input, levels, output;

        if (! cmd::loadInputVectorAndAdjustOthers(vm, input, output))
            return cmd::ERROR_UNHANDLED_EXCEPTION;
        BOOST_LOG_TRIVIAL(debug) << "Loaded input vector with " << input.size() << " samples.";

        if (! cmd::loadLevelsVector(vm, levels))
            return cmd::ERROR_UNHANDLED_EXCEPTION;
        BOOST_LOG_TRIVIAL(debug) << "Loaded levels vector with " << levels.size() << " elements.";
        
        VectorType data, weights;
        helpers::postprocessTVDNData(input, data, weights);
        BOOST_LOG_TRIVIAL(debug) << "Compressed input vector into " << data.size()
                                 << " (data, weight) tuples.";

        VectorType lambdas;
        cmd::loadLambdas(vm, lambdas);

        double prior_distance = vm.count("prior-distance") 
                                    ? vm["prior-distance"].as<double>()
                                    : 0.0;

        BinaryOptimization<> bin_opt(data.size(), levels.size());
        if (vm.count("debug-graphstructure"))
            bin_opt.recordEnergyGraphDumps();

        registerCostFunctions(bin_opt, 
                              data,
                              weights, 
                              levels, 
                              lambdas, 
                              prior_distance);

        auto energy = bin_opt.expansion(vm["maxiter"].as<int>());
        
        if (areAssignmentsRequested(vm)) {
            saveAssignments(vm, weights, levels, bin_opt);
            return cmd::SUCCESS;
        }

        postprocessAssignments(vm, bin_opt, input, weights, levels, output); 

        if (! cmd::saveOutputVector(vm, output))
            return cmd::ERROR_UNHANDLED_EXCEPTION;

        return cmd::SUCCESS;
    }
}

int main(const int argc, const char *argv[])
{
    using namespace std;

    auto desc = initializeOptionsDescription();

    try
    {
        // Define and parse the program options
        bpo::variables_map vm;
        if (! tryParseProgramOptions(desc, argc, argv, vm)) {
            return cmd::ERROR_IN_COMMAND_LINE;
        }

        cmd::configureLogging(vm.count("debug"));
        return runProgram(desc, vm);
    }
    catch(std::exception& e)
    {
        cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << endl;

        return cmd::ERROR_UNHANDLED_EXCEPTION;
    }

}
