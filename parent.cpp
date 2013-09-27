/*
 * File comments go here...
 */

#include "common.hpp"
#include <NHparser/NHparser.h>
#include <boost/program_options.hpp>
#include <iomanip>
#include <cstdlib>               // At least for EXIT_SUCCESS and EXIT_FAILURE
#include <utility>
#include <sstream>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>

//============================================================================
//             Namespace declarations and using directives.
//============================================================================

using namespace std;
using namespace boost::lambda;

//============================================================================
//                   Typedefs and class declarations
//============================================================================



//============================================================================
//                   Global Constants and variables.
//============================================================================

const string PROG_NAME = "parent";
const string USAGE = "Usage: " + PROG_NAME + " [OPTION]... FILE";

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
    typedef Rooted_tree< Param<double> > Hot;

    string filename;
    string orientation;

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    string opt;

    try
        {
            /* Declare visible options */
            visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("orientation,o",
                 po::value<string>(&orientation)->default_value("horizontal"),
                 "Possible values or 'horizontal' and 'vertical'")
                ;

            /* Declare positional options. */
            hidden_opts.add_options()
                ("filename",po::value<string>(&filename))
                ;

            po::positional_options_description positional_options;
            positional_options.add("filename", 1);

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
            cout
                << USAGE << '\n'
                << visible_opts << '\n';
            exit(EXIT_SUCCESS);
        }

    /* Check that all required positional arguments are given. */
    if (g_options.count("filename") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Check that orientation has been given correctly. */
    if (orientation != "horizontal" && orientation != "vertical")
        {
            cerr << PROG_NAME << ": Bad argument given to --orientation.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    FILE *infile = fopen(filename.c_str(), "r");
    if (errno)
        {
            cerr << PROG_NAME << ": could not open'" + filename + "' for reading\n";
            exit(EXIT_FAILURE);
        }

    NHNode *root = NH_read_tree(infile, stderr, 0);
    if (root == 0)
        {
            cerr << PROG_NAME << ": could not read tree from " + filename << "\n";
            exit(EXIT_FAILURE);
        }

    vector<NHNode *> nodes;
    nodes.push_back(root);
    for (unsigned i = 0; i < nodes.size(); ++i)
        {
            for (NHNodeList *l = nodes[i]->children; l != 0; l = l->next)
                {
                    nodes.push_back(l->node);
                }
        }

    BOOST_FOREACH (NHNode *n, nodes)
        {
            if (n->label == 0)
                {
                    cerr << PROG_NAME << ": some vertices lack label\n";
                    exit(EXIT_FAILURE);
                }
        }

    vector< pair<unsigned, unsigned> > edges;
    BOOST_FOREACH (NHNode *n, nodes)
        {
            unsigned u, v;
            istringstream in(n->label);
            in >> v;
            if (n->parent == 0)
                {
                    edges.push_back(make_pair(v, v));
                }
            else
                {
                    istringstream in(n->parent->label);
                    in >> u;
                    edges.push_back(make_pair(v, u));
                }
        }
    sort(edges.begin(), edges.end());

    unsigned width = 3;
    if (orientation == "horizontal")
        {
            for (unsigned i = 1; i < edges.size(); ++i)
                {
                    cout << setw(width) << edges[i].first;
                }
            cout << "\n";
            for (unsigned i = 1; i < edges.size(); ++i)
                {
                    cout << setw(width) << edges[i].second;
                }
            cout << "\n";
        }
    else
        {
            for (unsigned i = 1; i < edges.size(); ++i)
                {
                    cout << setw(width) << edges[i].first
                         << setw(width) << edges[i].second
                         << "\n";
                }
        }

    return EXIT_SUCCESS;
}

//=============================================================================
//                     Helper function declarations
//=============================================================================



//=============================================================================
//                   Function and member definitions.
//=============================================================================
