#ifndef COMMON_IMP_HPP
#define COMMON_IMP_HPP

#include "common.hpp"
#include <NHparser/NHparser.h>
#include <boost/foreach.hpp>
#include <sys/time.h>
#include <sstream>
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <stack>
#include <queue>
#include <list>
#include <set>
#include <map>

template <class RandomNumberGenerator>
void
init_rand(RandomNumberGenerator &engine)
{
    /*
     * Seed the random number generator.
     */
    timeval tv;
    gettimeofday(&tv, 0);
    unsigned seed = tv.tv_sec * 1000000 + tv.tv_usec;
    engine.seed(seed);
}

template <class RandomNumberGenerator>
void
init_rand(RandomNumberGenerator &engine, unsigned seed)
{
    /*
     * Seed the random number generator.
     */
    engine.seed(seed);
}

inline MP_double
log(MP_double d)
{
    return mpfr::Log(d);
}


template<typename Real>
void
gen_random_tree(unsigned N,
                unsigned root,
                std::vector<unsigned> &parent)
{
    /*
     * Get a random labeled tree by using the encoding of deo-02
     * (similar to prufer code but easier to decode)
     */

    /* Create a random encoding of length N-2. */
    std::vector<unsigned> code(N-2);
    for (unsigned i = 0; i < code.size(); ++i)
        {
            code[i] = g_rng_ui(N);
        }

    /* Create an adjacency list from the code. */
    std::stack<unsigned> S;
    std::vector<bool> used(N, false);

    for (int i = N-3; i >= 0; --i)
        {
            if (!used[code[i]])
                {
                    used[code[i]] = true;
                    S.push(code[i]);
                }
        }

    for (int i = N-1; i >=0; --i)
        {
            if (!used[i])
                S.push(i);
        }

    std::vector< std::list<unsigned> > neighbors(N);
    for (unsigned i = 0; i < N-2; ++i)
        {
            unsigned v = S.top(); S.pop();

            neighbors[code[i]].push_back(v);
            neighbors[v].push_back(code[i]);
        }

    unsigned a = S.top(); S.pop();
    unsigned b = S.top(); S.pop();

    neighbors[a].push_back(b);
    neighbors[b].push_back(a);

    /* Create the parent structure of the tree with 0 as root. */
    std::stack<unsigned> VS;
    VS.push(0);
    parent.resize(N);
    parent[0] = 0;

    while (!VS.empty())
        {
            unsigned u = VS.top(); VS.pop();
            BOOST_FOREACH (unsigned v, neighbors[u])
                {
                    if (v == parent[u])
                        continue;

                    parent[v] = u;
                    VS.push(v);
                }
        }
}


template<typename Real>
void
print_tree(std::ostream &out, Rooted_tree< Param<Real> > &tree, unsigned u)
{
    if (!tree.is_leaf(u))
        {
            out << "(";
            bool first = true;
            BOOST_FOREACH (unsigned child, tree.children(u))
                {
                    if (!first)
                        out << ", ";
                    print_tree(out, tree, child);
                    first = false;
                }
            out << ")";
        }
    out << u;
    out << "[&&"
        << tree.data(u).p_z_zp[0][0] << " "
        << tree.data(u).p_z_zp[0][1] << " "
        << tree.data(u).p_z_zp[1][0] << " "
        << tree.data(u).p_z_zp[1][1] << " "
        << tree.data(u).p_x_z[0][0] << " "
        << tree.data(u).p_x_z[0][1] << " "
        << tree.data(u).p_x_z[1][0] << " "
        << tree.data(u).p_x_z[1][1]
        << "]";
}


template<typename Real>
void
create_hot_from_NHNode(NHNode *root,
                       Rooted_tree< Param<Real> > &hot)
{
    using std::vector;
    using std::string;
    vector<NHNode *> nodes;

    /* Insert all nodes in the tree into the vector nodes. */
    nodes.push_back(root);
    for (unsigned i = 0; i < nodes.size(); ++i)
        {
            NHNode *n = nodes[i];
            for (NHNodeList *l = n->children; l != 0; l = l->next)
                {
                    nodes.push_back(l->node);
                }
        }

    /* Create the parent structure and check that labels are numeric
       and unique. Also extract the nodes' annotations. */
    vector<unsigned> parent(nodes.size());
    vector<bool> used_label(nodes.size());
    vector<string> annotation(nodes.size());
    unsigned root_number;
    BOOST_FOREACH (NHNode *n, nodes)
        {
            if (n->label == 0)
                throw Bad_hot("missing vertex label");

            /* Make sure label is numeric. */
            string s(n->label);
            if (s.find_first_not_of("0123456789") != s.npos)
                throw Bad_hot("bad vertex label (must be an integer)");

            /* Convert label to number. */
            std::istringstream in(s);
            unsigned u;
            in >> u;
            if (in.fail() || used_label[u] || u >= nodes.size())
                throw Bad_hot("bad vertex label");
            used_label[u] = true;

            if (n ->annotation == 0 || n->annotation->str == 0)
                throw Bad_hot("missing paramters");

            annotation[u] = n->annotation->str;

            unsigned p;
            if (n->parent == 0)
                {
                    root_number = u;
                    p = root_number;
                }
            else
                {
                    /* Convert parent's label to number. */
                    string s2 = n->parent->label;
                    std::istringstream in2(s2);
                    in2 >> p;
                }

            /* update the parent vector. */
            parent[u] = p;
        }

    hot.reset(parent, root_number);

    /* Set the parameters of the hot. */
    for (unsigned u = 0; u < parent.size(); ++u)
        {
            Param<Real> &param = hot.data(u);
            string s = annotation[u];
            s.erase(0, 3);
            s.erase(s.size() - 1, 1);

            std::istringstream in(s);
            in >> param.p_z_zp[0][0]
               >> param.p_z_zp[0][1]
               >> param.p_z_zp[1][0]
               >> param.p_z_zp[1][1]
               >> param.p_x_z[0][0]
               >> param.p_x_z[0][1]
               >> param.p_x_z[1][0]
               >> param.p_x_z[1][1];

            if (param.p_z_zp[0][0] < param.p_z_zp[1][0])
                param.p_z_zp[1][0] = Real(1) - param.p_z_zp[0][0];
            else
                param.p_z_zp[0][0] = Real(1) - param.p_z_zp[1][0];

            if (param.p_z_zp[0][1] < param.p_z_zp[1][1])
                param.p_z_zp[1][1] = Real(1) - param.p_z_zp[0][1];
            else
                param.p_z_zp[0][1] = Real(1) - param.p_z_zp[1][1];

            if (param.p_x_z[0][0] < param.p_x_z[1][0])
                param.p_x_z[1][0] = Real(1) - param.p_x_z[0][0];
            else
                param.p_x_z[0][0] = Real(1) - param.p_x_z[1][0];

            if (param.p_x_z[0][1] < param.p_x_z[1][1])
                param.p_x_z[1][1] = Real(1) - param.p_x_z[0][1];
            else
                param.p_x_z[0][1] = Real(1) - param.p_x_z[1][1];


            if (in.fail())
                throw Bad_hot("bad parameters");
        }
}


template<typename Real, typename Array_2d>
Real
compute_p_X(const Rooted_tree< Param<Real> > &hot,
            const Array_2d &data,
            unsigned dataset)
{
    boost::multi_array<Real, 2> p_Xu_u(boost::extents[hot.size()][2]);

    for (unsigned u = hot.postorder_begin(); u != hot.NONE; u = hot.postorder_next(u))
        {
            for (unsigned a = 0; a < 2; ++a)
                {
                    p_Xu_u[u][a] = hot.data(u).p_x_z[ data[dataset][u] ][a];

                    BOOST_FOREACH (unsigned v, hot.children(u))
                        {
                            p_Xu_u[u][a] *=
                                p_Xu_u[v][0] * hot.data(v).p_z_zp[0][a] +
                                p_Xu_u[v][1] * hot.data(v).p_z_zp[1][a];
                        }
                }
        }
//     for (unsigned u = 0; u < hot.size(); ++u)
//         std::cerr << data[dataset][u];
//     std::cerr << "\n";
//     for (unsigned u = 0; u < hot.size(); ++u)
//         {
//             std::cerr << p_Xu_u[u][0] << "\t" << p_Xu_u[u][1] << "\n";
//         }
//     for (unsigned u = 0; u < hot.size(); ++u)
//         {
//             for (unsigned a = 0; a < 2; ++a)
//                 {
//                     for (unsigned b = 0; b < 2; ++b)
//                         {
//                             std::cerr << hot.data(u).p_z_zp[a][b] << "\t";
//                         }
//                 }
//             for (unsigned a = 0; a < 2; ++a)
//                 {
//                     for (unsigned b = 0; b < 2; ++b)
//                         {
//                             std::cerr << hot.data(u).p_x_z[a][b] << "\t";
//                         }
//                 }
//             std::cerr << "\n";
//         }



    return
        hot.data(hot.root()).p_z_zp[0][1] * p_Xu_u[hot.root()][0] +
        hot.data(hot.root()).p_z_zp[1][1] * p_Xu_u[hot.root()][1];
}

template<typename Real>
void
read_mixture(std::string mixture_filename,
             std::vector< Rooted_tree< Param<Real> > > &hots,
             std::vector<Real> &probs)
{
    typedef Rooted_tree< Param<Real> > Hot;

    /* Open the file. */
    std::ifstream infile(mixture_filename.c_str());
    if (!infile)
        throw Bad_mixture("unable to open mixture file");

    /* Read the mixture file with filenames and probabilities. */
    std::vector<std::string> tree_files;
    probs.clear();

    std::string s;
    double d;
    infile >> s;
    if (infile.fail())
        {
            throw Bad_mixture("failed reading filename in mixture file");
        }
    while (infile)
        {
            infile >> d;
            if (infile.fail())
                {
                    throw Bad_mixture("failed reading probability in mixture file");
                }
            tree_files.push_back(s);
            probs.push_back(d);
            infile >> s;
        }

    /* Read the HOTs into the vector 'hots' */
    try
        {
            BOOST_FOREACH (std::string filename, tree_files)
                {
                    errno = 0;
                    FILE *infile = fopen(filename.c_str(), "r");
                    if (errno)
                        {
                            perror(0);
                            throw Bad_mixture("could not open '" + filename +
                                          "' for reading");
                        }

                    NHNode *root = NH_read_tree(infile, stderr, 0);
                    if (root == 0)
                        throw Bad_mixture("could not read HOT from '"
                                          + filename);

                    hots.push_back(Hot());
                    create_hot_from_NHNode(root, hots.back());
                    NH_destroy_tree(root);
                    fclose(infile);
                }
        }
    catch (Bad_hot b)
        {
            throw Bad_mixture(std::string("could not read HOT: ") + b.what());
        }

    /* Check that all HOTs have the same size and root */
    BOOST_FOREACH (Hot &hot, hots)
        {
            if (hot.size() != hots.back().size())
                throw Bad_mixture("hots in a mixture must all be of the same size");
            if (hot.root() != hots.back().root())
                throw Bad_mixture("hots in a mixture must all have the same vertex as root");
        }
}


// std::ostream &
// operator<<(std::ostream &out, const Log_double &ld)
// {
//     out << "'" << ld.raw_value() << "'";
//     if (ld == 0)
//         {
//             out << 0.0;
//         }
//     else
//         {
//             MP_double d = ld.raw_value();
//             d = mpfr::Exp(d);
//             out << d;
//         }
//     return out;
// }

#endif /* not COMMON_IMP_HPP */
