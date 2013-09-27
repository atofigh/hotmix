#include "common.hpp"
extern "C" {
#include <NHparser/NHparser.h>
}
#include "rooted-tree.hpp"
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <cstdlib>               // At least for EXIT_SUCCESS and EXIT_FAILURE
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include <utility>
#include <algorithm>
#include <sys/time.h>
#include <cstdio>
#include <cerrno>

//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace std;
namespace po = boost::program_options;

//============================================================================
//                   Typedefs and class declarations
//============================================================================

typedef Rooted_tree< Param<double> > Hot;

//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "gen-data";
const string USAGE = "Usage: " + PROG_NAME +
    " [OPTION]... MIXTURE_FILE NUM_DATASETS";

po::variables_map g_options;

//=============================================================================
//                        Function declarations
//=============================================================================

/*
 * gen_rand_bool()
 *
 * Returns true with probabilities p1 and zero with probability
 * p0. Uses the lesser of p0 and p1 to generate.
 */
bool gen_rand_bool(double p0, double p1)
{
    if (p0 < p1)
        return g_rng_d() >= p0;
    else
        return g_rng_d() < p1;
}


//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================

template<typename Real>
void
gen_data_from_hot(const Rooted_tree< Param<Real> > &hot,
                  vector<bool> &vec)
{
    vector<bool> Z(hot.size(), false);
    vec.resize(Z.size());

    for (unsigned u = hot.preorder_begin(); u != hot.NONE; u = hot.preorder_next(u))
        {
            if (u == hot.root())
                {
                    Z[u] = true;
                }
            else
                {
                    double r = g_rng_d();
                    Z[u] = r < hot.data(u).p_z_zp[1][ Z[hot.parent(u)] ];
                }
        }

    for (unsigned u = 0; u < Z.size(); ++u)
        {
            vec[u] = g_rng_d() < hot.data(u).p_x_z[1][ Z[u] ];
        }
}

//=============================================================================
//                                main()
//=============================================================================

int
main(int argc, char *argv[])
{
    string mixture_filename;
    unsigned num_datasets;
    bool visible_root = false;
    double p_root_z0 = -1;
    double p_root_x0_z0 = -1;
    double p_root_x0_z1 = -1;
    double p_root_x1_z0 = -1;
    double p_root_x1_z1 = -1;
    double p_z1_zp0 = -1;
    double p_x1_z0 = -1;

    init_rand(g_generator);

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    try
        {
            /* Declare visible options */
            visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("visible-root,r",
                 po::bool_switch(),
                 "If set data is also generated for the root")
                ("p-root-z0",
                 po::value<double>(&p_root_z0),
                 "Pr[Z(root)=0]")
                ("p-root-x0-z0",
                 po::value<double>(&p_root_x0_z0),
                 "Pr[X(root)=0|Z(root)=0]")
                ("p-root-x1-z0",
                 po::value<double>(&p_root_x1_z0),
                 "Pr[X(root)=1|Z(root)=0]")
                ("p-root-x0-z1",
                 po::value<double>(&p_root_x0_z1),
                 "Pr[X(root)=0|Z(root)=1]")
                ("p-root-x1-z1",
                 po::value<double>(&p_root_x1_z1),
                 "Pr[X(root)=1|Z(root)=1]")
                ("p-z1-zp0",
                 po::value<double>(&p_z1_zp0),
                 "Pr[Z(u)=1|Z(p(u))=0]")
                ("p-x1-z0",
                 po::value<double>(&p_x1_z0),
                 "Pr[X(u)=1|Z(u)=0]")
                ;

            /* Declare positional options. */
            hidden_opts.add_options()
                ("mixture-filename",po::value<string>(&mixture_filename))
                ("num-datasets",po::value<unsigned>(&num_datasets))
                ;

            po::positional_options_description positional_options;
            positional_options.add("mixture-filename", 1);
            positional_options.add("num-datasets", 1);

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
                 << visible_opts << '\n';
            exit(EXIT_SUCCESS);
        }

    /* Check that all required positional arguments are given. */
    if (g_options.count("num-datasets") == 0 ||
        g_options.count("mixture-filename") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }
    /* Check that the any probability given is in range. */
    if (g_options.count("p-root-z0") && (p_root_z0 < 0 || p_root_z0 > 1) ||
        g_options.count("p-root-x0-z0") && (p_root_x0_z0 < 0 || p_root_x0_z0 > 1) ||
        g_options.count("p-root-x0-z1") && (p_root_x0_z1 < 0 || p_root_x0_z1 > 1) ||
        g_options.count("p-root-x1-z0") && (p_root_x1_z0 < 0 || p_root_x1_z0 > 1) ||
        g_options.count("p-root-x1-z1") && (p_root_x1_z1 < 0 || p_root_x1_z1 > 1) ||
        g_options.count("p-z1-zp0") && (p_z1_zp0 < 0 || p_z1_zp0 > 1) ||
        g_options.count("p-x1-z0") && (p_x1_z0 < 0 || p_x1_z0 > 1))
        {
            cerr << PROG_NAME << ": Invalid probability given\n";
            exit(EXIT_FAILURE);
        }
    /* Check that at most one of any pair of complementary probs is given */
    if (g_options.count("p-root-x0-z0") && g_options.count("p-root-x1-z0") ||
        g_options.count("p-root-x0-z1") && g_options.count("p-root-x1-z1"))
        {
            cerr << PROG_NAME << ": At most one of a pair of complementary probabilities may be given\n";
            exit(EXIT_FAILURE);
        }
    /* Set the complementary probabilities. */
    if (g_options.count("p-root-x0-z0"))
        {
            p_root_x1_z0 = 1.0 - p_root_x0_z0;
        }
    if (g_options.count("p-root-x1-z0"))
        {
            p_root_x0_z0 = 1.0 - p_root_x1_z0;
        }
    if (g_options.count("p-root-x0-z1"))
        {
            p_root_x1_z1 = 1.0 - p_root_x0_z1;
        }
    if (g_options.count("p-root-x1-z1"))
        {
            p_root_x0_z1 = 1.0 - p_root_x1_z1;
        }


    /* Read the mixture. */
    vector<double> probs;
    vector<Hot> hots;
    try
        {
            read_mixture(mixture_filename, hots, probs);
        }
    catch (Bad_mixture b)
        {
            cerr << PROG_NAME << ": failed reading mixture file:\n";
            cerr << PROG_NAME << ": " << b.what() << "\n";
            exit(EXIT_FAILURE);
        }

    partial_sum(probs.begin(), probs.end(), probs.begin());

    /* check that the probabilities sum to ~1.0 */
    double prob_sum = probs.back();
    if (prob_sum > 1.0 + g_sum_tolerance || prob_sum < 1.0 - g_sum_tolerance)
        {
            cerr << PROG_NAME << ": mixture probabilites do not sum to ~1.0\n";
            exit(EXIT_FAILURE);
        }
    probs.back() = 1.0;

    /* Generate data from mixture. */
    for (unsigned i = 0; i < num_datasets; ++i)
        {
            /* choose a HOT from the mixture. */
            double r = g_rng_d();
            unsigned choice = 0;
            for (choice = 0; choice < probs.size(); ++choice)
                {
                    if (r < probs[choice])
                        break;
                }
            Hot &hot = hots[choice];
            unsigned root = hot.root();

            /* Generate from the chosen Hot. */
            vector<bool> vec(hot.size(), false);

            /* Generate the Z variables */
            vector<bool> Z(hot.size(), false);

            /* Generate Z for the root. */
            if (g_options.count("p-root-z0"))
                {
                    Z[root] = g_rng_d() >= p_root_z0;
                }
            else
                {
                    Z[root] = gen_rand_bool(hot.data(root).p_z_zp[0][1],
                                            hot.data(root).p_z_zp[1][1]);
                }

            /* Generate Z for the non-root vertices. */
            for (unsigned u = hot.preorder_begin(); u != hot.NONE; u = hot.preorder_next(u))
                {
                    if (u == root)
                        continue;

                    if (Z[hot.parent(u)])
                        {
                            Z[u] = gen_rand_bool(hot.data(u).p_z_zp[0][1],
                                                   hot.data(u).p_z_zp[1][1]);
                        }
                    else if (g_options.count("p-z1-zp0"))
                        {
                            Z[u] = g_rng_d() < p_z1_zp0;
                        }
                    else
                        {
                            Z[u] = gen_rand_bool(hot.data(u).p_z_zp[0][0],
                                                   hot.data(u).p_z_zp[1][0]);
                        }
                }

            /* Generate X for the root. */
            if (Z[root])
                {
                    if (g_options.count("p-root-x0-z1") ||
                        g_options.count("p-root-x1-z1"))
                        {
                            vec[root] = gen_rand_bool(p_root_x0_z1,
                                                      p_root_x1_z1);
                        }
                    else
                        {
                            vec[root]
                                = gen_rand_bool(hot.data(root).p_x_z[0][1],
                                                hot.data(root).p_x_z[1][1]);
                        }
                }
            else
                {
                    if (g_options.count("p-root-x0-z0") ||
                        g_options.count("p-root-x1-z0"))
                        {
                            vec[root] = gen_rand_bool(p_root_x0_z0,
                                                      p_root_x1_z0);
                        }
                    else
                        {
                            vec[root]
                                = gen_rand_bool(hot.data(root).p_x_z[0][0],
                                                hot.data(root).p_x_z[1][0]);
                        }
                }

            for (unsigned u = 0; u < hot.size(); ++u)
                {
                    if (u == root)
                        continue;

                    if (Z[u])
                        {
                            vec[u] = gen_rand_bool(hot.data(u).p_x_z[0][1],
                                                   hot.data(u).p_x_z[1][1]);
                        }
                    else if (g_options.count("p-x1-z0"))
                        {
                            vec[u] = g_rng_d() < p_x1_z0;
                        }
                    else
                        {
                            vec[u] = gen_rand_bool(hot.data(u).p_x_z[0][0],
                                                   hot.data(u).p_x_z[1][0]);
                        }
                }

            /* Print the generated dataset. */
            for (unsigned i = 0; i < root; ++i)
                {
                    cout << vec[i];
                }
            cout << " ";
            if (g_options["visible-root"].as<bool>())
                {
                    cout << vec[root] << " ";
                }
            for (unsigned i = root + 1; i < hot.size(); ++i)
                {
                    cout << vec[i];
                }
            cout << "\n";
        }


    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================



//=============================================================================
//                   Function and member definitions.
//=============================================================================
