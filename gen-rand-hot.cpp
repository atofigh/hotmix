#include "common.hpp"
#include "rooted-tree.hpp"
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <cstdlib>               // At least for EXIT_SUCCESS and EXIT_FAILURE
#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>

//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace std;
namespace po = boost::program_options;

//============================================================================
//                   Typedefs and class declarations
//============================================================================

typedef Rooted_tree< Param<double> > Hot;

struct Range {
public:
    struct Bad_range : public Exception {
        Bad_range(std::string message) : Exception(message) {}
    };

    Range(double m = 0, double M = 0);
    double m, M;
};
Range operator-(double d, Range &r);


//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "gen-rand-hot";
const string USAGE = "Usage: " + PROG_NAME + " [OPTION]... NUM_VERTICES ROOT";

po::variables_map g_options;

//=============================================================================
//                        Function declarations
//=============================================================================

void validate(boost::any& v,
              const std::vector<std::string>& values,
              Range* target_type, int);

//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================

//=============================================================================
//                                main()
//=============================================================================

int
main(int argc, char *argv[])
{
    unsigned num_vertices;
    unsigned root;
    int prec;
    Range default_same(0.5, 1.0);
    Range p_z_zp[2][2] = {{default_same, 1.0 - default_same},
                          {1.0 - default_same, default_same}};
    Range p_x_z[2][2] = {{default_same, 1.0 - default_same},
                          {1.0 - default_same, default_same}};

    Range default_root_same(0.9, 1.0);
    Range p_root_z[2] = {1.0 - default_root_same, default_root_same};
    Range p_root_x_z[2][2] = {{default_root_same, 1.0 - default_root_same},
                                {1.0 - default_root_same, default_root_same}};

    init_rand(g_generator);

    po::options_description phony_visible_opts("Command Line Options");
    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    try
        {
            /* Declare visible options */
            phony_visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("precision,p",
                 po::value<int>(&prec)->default_value(4),
                 "the number of digits with which the probabilities are written")
                ("p-zA-zpB",
                 "specifies the range for the probability\nPr[Z(u)=A|Z(p(u))=B], "
                 "where A and B are either 0 or 1.")
                ("p-xA-zB",
                 "specifies the range for the probability\nPr[X(u)=A|Z(u)=B], "
                 "where A and B are either 0 or 1.")
                ("p-root-zA",
                 "specifies the range for the probability\nPr[Z(root)=A], "
                 "where A is either 0 or 1.")
                ("p-root-xA_zB",
                 "specifies the range for the probability\nPr[X(root)=A|Z(root)=B], "
                 "where A and B are either 0 or 1.")
                ("common-error-probs",
                 "If set, all edges will share a single probability for"
                 "Pr[X(u)=1|Z(u)=0], Pr[X(u)=0|Z(u)=1], and Pr[Z(u)=1|Z(p(u))=0].")
                ;
            visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("precision,p",
                 po::value<int>(&prec)->default_value(4),
                 "the number of digits with which the probabilities are written")
                ("p-z0-zp0", po::value<Range>(&p_z_zp[0][0]), "")
                ("p-z0-zp1", po::value<Range>(&p_z_zp[0][1]), "")
                ("p-z1-zp0", po::value<Range>(&p_z_zp[1][0]), "")
                ("p-z1-zp1", po::value<Range>(&p_z_zp[1][1]), "")
                ("p-x0-z0", po::value<Range>(&p_x_z[0][0]), "")
                ("p-x0-z1", po::value<Range>(&p_x_z[0][1]), "")
                ("p-x1-z0", po::value<Range>(&p_x_z[1][0]), "")
                ("p-x1-z1", po::value<Range>(&p_x_z[1][1]), "")
                ("p-root-z0", po::value<Range>(&p_root_z[0]), "")
                ("p-root-z1", po::value<Range>(&p_root_z[1]), "")
                ("p-root-x0-z0", po::value<Range>(&p_root_x_z[0][0]), "")
                ("p-root-x0-z1", po::value<Range>(&p_root_x_z[0][1]), "")
                ("p-root-x1-z0", po::value<Range>(&p_root_x_z[1][0]), "")
                ("p-root-x1-z1", po::value<Range>(&p_root_x_z[1][1]), "")
                ("common-error-probs",
                 po::bool_switch(),
                 "If set, all edges will share a single probability for"
                 "Pr[X(u)=1|Z(u)=0], Pr[X(u)=0|Z(u)=1], and Pr[Z(u)=1|Z(p(u))=0].")
                ;

            /* Declare positional options. */
            hidden_opts.add_options()
                ("num-vertices", po::value<unsigned>(&num_vertices))
                ("root", po::value<unsigned>(&root))
                ;

            po::positional_options_description positional_options;
            positional_options.add("num-vertices", 1);
            positional_options.add("root", 1);

            /* Gather all options into a single options_description. */
            all_options.add(visible_opts).add(hidden_opts);

            /* Parse the arguments. */
            po::command_line_parser parser(argc, argv);
            parser.options(all_options);
            parser.positional(positional_options);
            po::store(parser.run(), g_options);
            po::notify(g_options);
        }
    catch (exception &e)
        {
            cerr << PROG_NAME << ": " << e.what() << '\n'
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Show help message if --help was given. */
    if (g_options.count("help"))
        {
            cout << USAGE << '\n'
                 << phony_visible_opts << '\n';
            exit(EXIT_SUCCESS);
        }

    /* Check that all required positional arguments are given. */
    if (g_options.count("num-vertices") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }
    /* Check that num_vertices > 1. */
    if (num_vertices <= 1)
        {
            cerr << PROG_NAME << ": number of vertices must be at least 2.\n";
            exit(EXIT_FAILURE);
        }
    /* Check that root < num_vertices. */
    if (root >= num_vertices)
        {
            cerr << PROG_NAME << ": ROOT must be less than NUM_VERTICES.\n";
            exit(EXIT_FAILURE);
        }
    /* Check that at most on of p-#0-#A and p-#1-#A is given. */
    if (g_options.count("p-z0-zp0") + g_options.count("p-z1-zp0") > 1 ||
        g_options.count("p-z0-zp1") + g_options.count("p-z1-zp1") > 1 ||
        g_options.count("p-x0-z0") + g_options.count("p-x1-z0") > 1 ||
        g_options.count("p-x0-z1") + g_options.count("p-x1-z1") > 1 ||
        g_options.count("p-root-z0") + g_options.count("p-root-z1") > 1 ||
        g_options.count("p-root-x0-z0)") + g_options.count("p-root-x1-z0)") > 1 ||
        g_options.count("p-root-x0-z1)") + g_options.count("p-root-x1-z1)") > 1)
        {
            cerr << PROG_NAME << ": incompatible probability ranges given.\n";
            exit(EXIT_FAILURE);
        }
    /* set the ranges. */
    if (g_options.count("p-z0-zp0"))
        {
            p_z_zp[0][0] = g_options["p-z0-zp0"].as<Range>();
            p_z_zp[1][0] = 1.0 - p_z_zp[0][0];
        }
    if (g_options.count("p-z1-zp0"))
        {
            p_z_zp[1][0] = g_options["p-z1-zp0"].as<Range>();
            p_z_zp[0][0] = 1.0 - p_z_zp[1][0];
        }
    if (g_options.count("p-z0-zp1"))
        {
            p_z_zp[0][1] = g_options["p-z0-zp1"].as<Range>();
            p_z_zp[1][1] = 1.0 - p_z_zp[0][1];
        }
    if (g_options.count("p-z1-zp1"))
        {
            p_z_zp[1][1] = g_options["p-z1-zp1"].as<Range>();
            p_z_zp[0][1] = 1.0 - p_z_zp[1][1];
        }
    if (g_options.count("p-x0-z0"))
        {
            p_x_z[0][0] = g_options["p-x0-z0"].as<Range>();
            p_x_z[1][0] = 1.0 - p_x_z[0][0];
        }
    if (g_options.count("p-x1-z0"))
        {
            p_x_z[1][0] = g_options["p-x1-z0"].as<Range>();
            p_x_z[0][0] = 1.0 - p_x_z[1][0];
        }
    if (g_options.count("p-x0-z1"))
        {
            p_x_z[0][1] = g_options["p-x0-z1"].as<Range>();
            p_x_z[1][1] = 1.0 - p_x_z[0][1];
        }
    if (g_options.count("p-x1-z1"))
        {
            p_x_z[1][1] = g_options["p-x1-z1"].as<Range>();
            p_x_z[0][1] = 1.0 - p_x_z[1][1];
        }
    if (g_options.count("p-root-z0"))
        {
            p_root_z[0] = g_options["p-root-z0"].as<Range>();
            p_root_z[1] = 1.0 - p_root_z[0];
        }
    if (g_options.count("p-root-z1"))
        {
            p_root_z[1] = g_options["p-root-z1"].as<Range>();
            p_root_z[0] = 1.0 - p_root_z[1];
        }
    if (g_options.count("p-root-x0-z0"))
        {
            p_root_x_z[0][0] = g_options["p-root-x0-z0"].as<Range>();
            p_root_x_z[1][0] = 1.0 - p_root_x_z[0][0];
        }
    if (g_options.count("p-root-x1-z0"))
        {
            p_root_x_z[1][0] = g_options["p-root-x1-z0"].as<Range>();
            p_root_x_z[0][0] = 1.0 - p_root_x_z[1][0];
        }
    if (g_options.count("p-root-x0-z1"))
        {
            p_root_x_z[0][1] = g_options["p-root-x0-z1"].as<Range>();
            p_root_x_z[1][1] = 1.0 - p_root_x_z[0][1];
        }
    if (g_options.count("p-root-x1-z1"))
        {
            p_root_x_z[1][1] = g_options["p-root-x1-z1"].as<Range>();
            p_root_x_z[0][1] = 1.0 - p_root_x_z[1][1];
        }

    /* Get a random tree. */
    Hot hot;
    vector<unsigned> parent;
    gen_random_tree<double>(num_vertices, root, parent);
    hot.reset(parent, root);

    /* Set the parameters according to the ranges. */

    Param<double> common_param;
    for (unsigned b = 0; b < 2; ++b)
        {
            unsigned a = p_z_zp[0][b].m < p_z_zp[1][b].m ? 0 : 1;
            double m = p_z_zp[a][b].m;
            double M = p_z_zp[a][b].M;

            common_param.p_z_zp[a][b] = g_rng_d() * (M - m) + m;
            common_param.p_z_zp[1-a][b] = 1.0 - common_param.p_z_zp[a][b];
        }
    for (unsigned b = 0; b < 2; ++b)
        {
            unsigned a = p_x_z[0][b].m < p_x_z[1][b].m ? 0 : 1;
            double m = p_x_z[a][b].m;
            double M = p_x_z[a][b].M;

            common_param.p_x_z[a][b] = g_rng_d() * (M - m) + m;
            common_param.p_x_z[1-a][b] = 1.0 - common_param.p_x_z[a][b];
        }


    for (unsigned u = 0; u < hot.size(); ++u)
        {
            Param<double> &param = hot.data(u);

            for (unsigned b = 0; b < 2; ++b)
                {
                    unsigned a = p_z_zp[0][b].m < p_z_zp[1][b].m ? 0 : 1;
                    double m = p_z_zp[a][b].m;
                    double M = p_z_zp[a][b].M;

                    param.p_z_zp[a][b] = g_rng_d() * (M - m) + m;
                    param.p_z_zp[1-a][b] = 1.0 - param.p_z_zp[a][b];
                }
            for (unsigned b = 0; b < 2; ++b)
                {
                    unsigned a = p_x_z[0][b].m < p_x_z[1][b].m ? 0 : 1;
                    double m = p_x_z[a][b].m;
                    double M = p_x_z[a][b].M;

                    param.p_x_z[a][b] = g_rng_d() * (M - m) + m;
                    param.p_x_z[1-a][b] = 1.0 - param.p_x_z[a][b];
                }
            if (g_options["common-error-probs"].as<bool>())
                {
                    param.p_z_zp[0][0] = common_param.p_z_zp[0][0];
                    param.p_z_zp[1][0] = common_param.p_z_zp[1][0];
                    param.p_x_z[0][0] = common_param.p_x_z[0][0];
                    param.p_x_z[0][1] = common_param.p_x_z[0][1];
                    param.p_x_z[1][0] = common_param.p_x_z[1][0];
                    param.p_x_z[1][1] = common_param.p_x_z[1][1];
                }
        }
    Param<double> &param = hot.data(hot.root());

    unsigned a = p_root_z[0].m < p_root_z[1].m ? 0 : 1;
    double m = p_root_z[a].m;
    double M = p_root_z[a].M;

    param.p_z_zp[a][1] = g_rng_d() * (M - m) + m;
    param.p_z_zp[1-a][1] = 1.0 - param.p_z_zp[a][1];
    param.p_z_zp[0][0] = 0;
    param.p_z_zp[1][0] = 0;
    for (unsigned b = 0; b < 2; ++b)
        {
            unsigned a = p_root_x_z[0][b].m < p_root_x_z[1][b].m ? 0 : 1;
            double m = p_root_x_z[a][b].m;
            double M = p_root_x_z[a][b].M;

            param.p_x_z[a][b] = g_rng_d() * (M - m) + m;
            param.p_x_z[1-a][b] = 1.0 - param.p_x_z[a][b];
        }

    cout.precision(prec);
    print_tree(cout, hot, hot.root());
    cout << ";\n";

    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================



//=============================================================================
//                   Function and member definitions.
//=============================================================================

void validate(boost::any& v,
              const std::vector<std::string>& values,
              Range* target_type, int)
{
    po::validators::check_first_occurrence(v);
    const string& s = po::validators::get_single_string(values);
    double m, M;

    istringstream in(s);
    char c;
    string k;
    in >> m;
    in >> c;
    in >> M
       >> ws;

    if (in.fail() || c != ',' || !in.eof())
        throw po::invalid_option_value("invalid value");

    try
        {
            v = boost::any(Range(m, M));
        }
    catch (Range::Bad_range b)
        {
            throw po::invalid_option_value(b.what());
        }
}

Range::Range(double m, double M) : m(m), M(M)
{
    if (m > M ||
        m < 0 || m > 1 ||
        M < 0 || M > 1)
        {
            throw Bad_range("invalid range");
        }
}

Range
operator-(double d, Range &r)
{
    double M = d - r.m;
    double m = d - r.M;
    return Range(m, M);
}
