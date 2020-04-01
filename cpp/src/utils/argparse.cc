#include <boost/program_options.hpp>

// TODO
namespace prog_opts = boost::program_options;

// prog_opts::options_description cli_opts;

/*
    po::options_description gen_opts("General options");

    gen_opts.add_options()
        ("help,h", "show help message")
        ("pool", po::value<std::string>(&pool)->required(), "pool")
    ;

    po::options_description all_opts("Allowed options");
    all_opts.add(gen_opts);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opts), vm);

    if (vm.count("help")) {
        std::cout << all_opts << std::endl;
        return 1;
    }

    po::notify(vm);
*/
