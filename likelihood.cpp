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

//============================================================================
//                   Typedefs and class declarations
//============================================================================

//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "likelihood";
const string USAGE = "Usage: " + PROG_NAME +
    " [OPTION]... MIXTURE_FILE DATA_FILE";

boost::program_options::variables_map g_options;

//=============================================================================
//                        Function declarations
//=============================================================================



//=============================================================================
//         Template and inline function and member definitions.
//=============================================================================



//=============================================================================
//                                main()
//=============================================================================

int
main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    typedef Rooted_tree< Param<MP_double> > Hot;


    string mixture_filename;
    string data_filename;
    double p_root_z0 = -1;
    double p_root_x0_z0 = -1;
    double p_root_x0_z1 = -1;
    double p_root_x1_z0 = -1;
    double p_root_x1_z1 = -1;
    double p_z1_zp0 = -1;
    double p_x1_z0 = -1;

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    try
        {
            /* Declare visible options */
            visible_opts.add_options()
                ("help,h", "display this help and exit")
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
                ("data-filename",po::value<string>(&data_filename))
                ;

            po::positional_options_description positional_options;
            positional_options.add("mixture-filename", 1);
            positional_options.add("data-filename", 1);

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
    if (g_options.count("data-filename") == 0 ||
        g_options.count("mixture-filename") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Check that p_z1_zp0 and p_x1_z0 are reasonable. */
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
    vector<MP_double> probs;
    vector<Hot> hots;
    try
        {
            read_mixture<MP_double>(mixture_filename, hots, probs);
        }
    catch (Bad_mixture b)
        {
            cerr << PROG_NAME << ": failed reading mixture file '" << mixture_filename << "':\n";
            cerr << PROG_NAME << ": " << b.what() << "\n";
            exit(EXIT_FAILURE);
        }

    /* Set the special probabilities of the hots as given on cmd-line. */
    const unsigned root = hots.back().root();
    BOOST_FOREACH (Hot &hot, hots)
        {
            if (p_root_z0 >= 0)
                {
                    hot.data(root).p_z_zp[0][1] = p_root_z0;
                    hot.data(root).p_z_zp[1][1] = 1.0 - p_root_z0;
                }
            if (p_root_x0_z0 >= 0)
                {
                    hot.data(root).p_x_z[0][0] = p_root_x0_z0;
                    hot.data(root).p_x_z[1][0] = p_root_x1_z0;
                }
            if (p_root_x0_z1 >= 0)
                {
                    hot.data(root).p_x_z[0][1] = p_root_x0_z1;
                    hot.data(root).p_x_z[1][1] = p_root_x1_z1;
                }
            if (p_z1_zp0 >= 0)
                {
                    for (unsigned u = 0; u != hot.size(); ++u)
                        {
                            if (u == root)
                                continue;

                            hot.data(u).p_z_zp[1][0] = p_z1_zp0;
                            hot.data(u).p_z_zp[0][0] = 1.0 - p_z1_zp0;
                        }
                }
            if (p_x1_z0 >= 0)
                {
                    for (unsigned u = 0; u != hot.size(); ++u)
                        {
                            if (u == root)
                                continue;

                            hot.data(u).p_x_z[1][0] = p_x1_z0;
                            hot.data(u).p_x_z[0][0] = 1.0 - p_x1_z0;
                        }
                }
        }

    /* check that the mixture probabilities sum to ~1.0 */
    MP_double prob_sum = accumulate(probs.begin(), probs.end(), MP_double(0.0));
    if (prob_sum > MP_double(1.0 + g_sum_tolerance) ||
        prob_sum < MP_double(1.0 - g_sum_tolerance))
        {
            cerr << PROG_NAME << ": mixture probabilites do not sum to ~1.0\n";
            exit(EXIT_FAILURE);
        }

    /* Read the data. */
    ifstream infile(data_filename.c_str());
    if (!infile)
        {
            cerr << PROG_NAME << ": Could not open file '" << data_filename<< "'\n";
            exit(EXIT_FAILURE);
        }

    vector<string> input_lines;
    read_input_lines(infile, input_lines);
    infile.close();

    /* Do we have at least one input line? */
    if (input_lines.empty())
        {
            cerr << PROG_NAME << ": no valid input lines detected.\n";
            exit(EXIT_FAILURE);
        }

    unsigned num_datasets = input_lines.size();
    unsigned num_columns = input_lines.front().size();

    /* Do all lines have the same number of variables? */
    BOOST_FOREACH (string s, input_lines)
        {
            if (s.size() != num_columns)
                {
                    cerr << PROG_NAME << ": error reading input file\n";
                    cerr << PROG_NAME << ": every valid input line must contain the same number of variables\n";
                    exit(EXIT_FAILURE);
                }
        }

    /* Is the number of variables equal to the size of the trees? */
    if (num_columns != hots.back().size() &&
        num_columns != hots.back().size() - 1)
        {
            cerr << PROG_NAME << ": error reading mixture\n";
            cerr << PROG_NAME << ": the number of data columns does not match the size of the trees\n";
            exit(EXIT_FAILURE);
        }

    /* Create the data matrix. */
    boost::multi_array<bool, 2> data(boost::extents[num_datasets][hots.back().size()]);
    unsigned r = 0;
    BOOST_FOREACH (string s, input_lines)
        {
            unsigned c = 0;
            BOOST_FOREACH (char ch, s)
                {
                    if (c == root && num_columns < hots.back().size())
                        {
                            data[r][c] = true;
                            c++;
                        }
                    data[r][c] = ch == '0' ? false : true;
                    ++c;
                }
            ++r;
        }
    clog << "Read " << num_datasets << " datasets with "
         << num_columns << " columns in each.\n";

    /* compute the likelihood of the mixture. */
    MP_double L(1.0);
    for (unsigned d = 0; d < num_datasets; ++d)
        {
            MP_double factor(0.0);
            for (unsigned k = 0; k < hots.size(); ++k)
                {
                    factor +=
                        probs[k] * compute_p_X<MP_double>(hots[k], data, d);
                }
            L *= factor;
        }

    cout << L << "\n";

    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================



//=============================================================================
//                   Function and member definitions.
//=============================================================================
