#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/log/trivial.hpp>

#include "common/cmd_helpers.h"


#include "helpers.h"
#include "condat_denoise.h"

namespace bpo = boost::program_options;
namespace ublas = boost::numeric::ublas;

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
                        "containing the noisiy input data set")
                ("output", bpo::value<string>()->default_value("-"),
                        "Filename of the matrix market vector file "
                        "the denoised data should be written to")
                ("lambda", bpo::value<double>(),
                        "Lambda coefficient used as regularizer in"
                        "the total-variation denoising problem")
                ("debug,d", "Turn on debug output if flag is set");

        return desc;
    }

    bool tryParseProgramOptions(const bpo::options_description& desc,
                                const int argc,
                                const char *argv[],
                                bpo::variables_map& vm)
    {
        using namespace std;
        try {
            bpo::positional_options_description ppos;
            ppos.add("input", 1);
            ppos.add("output", 2);

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
                BOOST_LOG_TRIVIAL(error) << "'input' argument is required";
                is_valid = false;
            }

            if (! vm.count("output")) {
                BOOST_LOG_TRIVIAL(error) << "'output' argument is required";
                is_valid = false;
            }

            if (! vm.count("lambda")) {
                BOOST_LOG_TRIVIAL(error) << "'lambda' argument is required";
                is_valid = false;
            }

            return is_valid;
        }
        catch (bpo::error& e) {
            BOOST_LOG_TRIVIAL(error) << e.what();
            BOOST_LOG_TRIVIAL(error) << desc;
            return false;
        }
    }

    int runProgram(const bpo::options_description& desc,
                   const bpo::variables_map &vm)
    {
        using namespace std;
        typedef boost::numeric::ublas::vector<double> VectorType;

	auto lambda = vm["lambda"].as<double>();

        VectorType input, output;
        if (! cmd::loadInputVectorAndAdjustOutputSize(vm, input, output))
            return cmd::ERROR_UNHANDLED_EXCEPTION;

        TV1D_denoise(input, output, lambda);

        if (! cmd::saveOutputVector(vm, output))
            return cmd::ERROR_UNHANDLED_EXCEPTION;

        return cmd::SUCCESS;
    }
}


int main(const int argc, const char *argv[])
{
    using namespace std;

    cmd::configureLogging();
    auto desc = initializeOptionsDescription();

    try
    {
        // Define and parse the program options
        bpo::variables_map vm;
        if (! tryParseProgramOptions(desc, argc, argv, vm)) {
            return cmd::ERROR_IN_COMMAND_LINE;
        }

        return runProgram(desc, vm);
    }
    catch(std::exception& e)
    {

        cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << endl;

        return cmd::ERROR_UNHANDLED_EXCEPTION;
    }
}
