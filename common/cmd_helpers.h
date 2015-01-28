#ifndef CMD_HELPER_H
#define CMD_HELPER_H

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/expressions.hpp>
#include <boost/core/null_deleter.hpp>

#include "helpers.h"

#ifndef NAN
    static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
    #define NAN (*(const float *) __nan)
#endif
	

namespace bpo = boost::program_options;

namespace cmd
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

	
    void configureLogging(bool debug = false)
    {
        using namespace std;
        using namespace boost::log;

        typedef sinks::synchronous_sink<sinks::text_ostream_backend> text_sink;
        boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
        boost::shared_ptr<std::ostream> stream(&std::clog, boost::null_deleter());

        sink->locked_backend()->add_stream(stream);
        sink->set_formatter
                (
                        expressions::format("<%1%> %2%")
                                % trivial::severity
                                % expressions::smessage
                );
	if (! debug) {
	    sink->set_filter(trivial::severity >= trivial::info);
	}

        core::get()->remove_all_sinks();
        core::get()->add_sink(sink);
    }

    bool isHelpRequest(const bpo::options_description &desc,
                       const bpo::variables_map &vm)
    {
        using namespace std;

        // --help option
        if (vm.count("help"))
        {
            cout << "graph_processing command line interface" << endl
                    << desc << endl;

            return true;
        }

        return false;
    }

    template<typename VectorType>
    bool loadVector(const std::string &filename,
		    VectorType& v)
    {
	using namespace std;

	if (! helpers::loadMMVector(v, filename)) {
     	    BOOST_LOG_TRIVIAL(error) << "Error during loading of vector "
                        << "from file '" << filename << "'";
	    return false;
	}

	return true;
    }

    template <typename VectorType>
    bool loadInputVector(const bpo::variables_map &vm,
                         VectorType& input)
    {
        using namespace std;

        return loadVector(vm["input"].as<string>(), input);
    }

    template <typename VectorType>
    bool loadInputVectorAndAdjustOutputSize(const bpo::variables_map &vm,
                                            VectorType& input,
                                            VectorType& output)
    {
        if (! loadInputVector(vm, input))
            return false;

        output.resize(input.size());
        return true;
    }

	template <typename VectorType>
    bool loadInputVectorAndAdjustOthers(const bpo::variables_map &vm,
		    			VectorType &input,
						VectorType &o1)
    {
		if (! loadInputVector(vm, input))
			return false;
		
		o1.resize(input.size());
		
		return true;
    }
	
    template <typename VectorType>
    bool loadInputVectorAndAdjustOthers(const bpo::variables_map &vm,
		    			VectorType &input,
						VectorType &o1,
						VectorType &o2)
    {
		if (! loadInputVector(vm, input))
			return false;
		
		o1.resize(input.size());
		o2.resize(input.size());
		
		return true;
    }
	
	
    template <typename VectorType>
    bool loadLevelsVector(const bpo::variables_map &vm,
		          VectorType &levels)
    {
		using namespace std;

		return loadVector(vm["levels"].as<string>(), levels);
    }

    template <typename VectorType>
    void loadLambdas(const bpo::variables_map &vm,
                     VectorType &lambdas)
    {
        double rho_d = vm["rho-d"].as<double>();
        double rho_s = vm["rho-s"].as<double>();
        double rho_p = vm["rho-p"].as<double>();

        lambdas.resize(3);
        lambdas(0) = rho_d != 0.0 ? (1.0/rho_d) : 0.0;
        lambdas(1) = rho_s != 0.0 ? (1.0/rho_s) : 0.0;
        lambdas(2) = rho_p != 0.0 ? (1.0/rho_p) : 0.0;
    }

    template <typename VectorType>
    void loadJumpDistParams(const bpo::variables_map &vm,
		    	    VectorType &jump_dist_params)
    {
	jump_dist_params.resize(3);
	jump_dist_params(0) = 0.0;
	jump_dist_params(1) = 0.0;
	jump_dist_params(2) = 0.0;
    }

    std::shared_ptr<std::ostream> openOutputStream(const bpo::variables_map &vm)
    {
        using namespace std;

        if(vm["output"].as<string>() == "-") {
            return shared_ptr<ostream>(&cout, [](ostream *) {});
        }

        return make_shared<ofstream>(vm["output"].as<string>());
    }

    template <typename VectorType>
    bool saveOutputVector(const bpo::variables_map &vm,
                          const VectorType &output)
    {
        using namespace std;

        auto os = openOutputStream(vm);
        return helpers::saveMMVector(output, *os);
    }

    template <typename MatrixType>
    bool saveOutputMatrix(const bpo::variables_map &vm,
		      	  const MatrixType &output)
    {
	using namespace std;

	auto os = openOutputStream(vm);
	return helpers::saveMMMatrix(output, *os);
    }

    double convertToDouble(const std::string& s)
    {
		using namespace std;
		
        try
        {
            return boost::lexical_cast<double>(s);
        }
        catch (boost::bad_lexical_cast const&)
        {
            return NAN;
        }
    }

    inline bool tryLoadDataAndWeightsLine(const std::string &line,
            std::pair<double, double>& p)
    {
        using namespace std;

        boost::char_separator<char> seps(", ");
        boost::tokenizer<boost::char_separator<char>> tokenizer(line, seps);
        auto it = tokenizer.begin();

        auto t1 = *it++;
        auto t2 = *it;

        boost::algorithm::trim(t1);
        boost::algorithm::trim(t2);

        p = make_pair(convertToDouble(t1),
                convertToDouble(t2));

        return true;
    }

    bool tryLoadDataAndWeights(const std::string& filename,
            std::vector<double>& data,
            std::vector<double>& weights)
    {
        using namespace std;

        data.clear();
        weights.clear();

        ifstream in_file(filename);
        if (! in_file.is_open()) {
            cerr << "ERROR: Unable to open file '" << filename << "'" << endl;
            return false;
        }

        string current_line;
        pair<double, double> current_pair;
        while (getline(in_file, current_line))
        {
            if (! tryLoadDataAndWeightsLine(current_line, current_pair))
                continue;

            data.push_back(get<0>(current_pair));
            weights.push_back(get<1>(current_pair));
        }

        if (data.size() == 0 || weights.size() == 0) {
            cerr << "ERROR: The file '" << filename << "' was empty" << endl;
            return false;
        }

        return true;
    }

    bool tryLoadLabels(const std::string& filename,
            std::vector<double>& levels)

    {
        using namespace std;

        levels.clear();

        ifstream in_file(filename);
        if (! in_file.is_open()) {
            cerr << "ERROR: Unable to open file '" << filename << "'" << endl;
            return false;
        }

        string current_line;
        while (getline(in_file, current_line))
        {
            boost::algorithm::trim(current_line);
            levels.push_back(convertToDouble(current_line));
        }

        if (levels.size() == 0) {
            cerr << "ERROR: The file '" << filename << "' was empty" << endl;
            return false;
        }

        return true;
    }

    bool tryLoadLambdas(const std::string& filename,
            std::vector<double>& levels)
    {
        return tryLoadLabels(filename, levels);
    }


}

#endif
