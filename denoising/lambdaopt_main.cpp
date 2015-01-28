#include <limits>
#include <memory>

#include "common/cmd_helpers.h"
#include "lambda_max.h"
#include "lambda_opt.h"

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
                ("lambdamax", "Just output lambda max, "
                              "which is the maximum value of the regulrization parameter")
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

            return is_valid;
        }
        catch (bpo::error& e) {
            BOOST_LOG_TRIVIAL(error) << e.what();
            BOOST_LOG_TRIVIAL(error) << desc;
            return false;
        }
    }

    bool isLambdaMaxRequest(const bpo::options_description &desc,
                            const bpo::variables_map &vm)
    {
        using namespace std;

        if (vm.count("lambdamax"))
            return true;

        return false;
    }

    int runProgram(const bpo::options_description& desc,
                   const bpo::variables_map &vm)
    {
        using namespace std;
        typedef boost::numeric::ublas::vector<double> VectorType;

        VectorType input;
        if (! cmd::loadInputVector(vm, input))
            return cmd::ERROR_UNHANDLED_EXCEPTION;

	double lambda_min = 100; 
        double lambda_max = computeLambdaMax(input);
        if (isLambdaMaxRequest(desc, vm)) {
            cout << lambda_max << endl;
            return cmd::SUCCESS;
        }


        double lambda_opt = computeLambdaOpt(input, lambda_min, lambda_max);
        cout << lambda_opt << endl;


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
