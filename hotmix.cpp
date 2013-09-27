/*
 * File comments go here...
 */

#include "common.hpp"
#include "edmonds_optimum_branching.hpp"
#include "rooted-tree.hpp"
#include <mpfrcpp/mpfrcpp.hpp>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/operators.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <limits>
#include <sys/time.h>
#include <cstdio>
#include <cerrno>
extern "C" {
#include <NHparser/NHparser.h>
}

//*****************************************************************************
//             Namespace declarations and using directives.
//*****************************************************************************

namespace        po = boost::program_options;
using namespace  std;
using            boost::multi_array;
using            boost::dynamic_bitset;

//*****************************************************************************
//                   Typedefs and class declarations
//*****************************************************************************

typedef MP_double                           real_type;
typedef Rooted_tree< Param<real_type> >     Hot;

/*
 * Class Program_input
 *
 * A global object 'g_in' of Program_input is used to hold input to
 * the program, and can be made accessible to any functions/modules
 * that would benefit from it.
 */
struct Program_input {
    multi_array<bool, 2>        data;
    int                         num_datasets;
    int                         num_variables;
    vector<unsigned>            input_freq;
    unsigned                    total_datasets; // including copies

    unsigned                    num_trees;
    vector<string>              start_trees;
    vector<double>              start_lambdas;

    vector< dynamic_bitset<> >  column_masks;

    int                         max_mixture_iterations;
    int                         max_tree_iterations;
    double                      mixture_cutoff;
    double                      tree_cutoff;

    bool                        const_tree;
    bool                        global_x_params;

    double                      pseudo_count;
    double                      p_z1_zp0;
    double                      p_x1_z0;
};

//*****************************************************************************
//                   Global Constants and variables.
//*****************************************************************************

const string PROG_NAME = "hotmix";
const string USAGE = "Usage: " + PROG_NAME + " [OPTION]... DATA_FILE NUM_TREES";

/* Variable holding the command line options. */
po::variables_map g_options;

/* The global variable holding the program input. */
Program_input g_in;

/* Variables holding the hots and gammas */
vector< vector<real_type> > g_gamma;
vector<Hot>                 g_trees;
vector<real_type>           g_lambda;
unsigned                    g_root = 0;

/* The definition of a complete graph for use with BGL-style
   edmonds-alg. */
struct complete_graph {
    int n_vertices;

    complete_graph(int n_vertices) : n_vertices(n_vertices) {}

    struct edge_iterator : public boost::input_iterator_helper<edge_iterator, int, void, void, int>
    {
        int edge_idx, n_vertices;

        edge_iterator() : edge_idx(0), n_vertices(-1) {}
        edge_iterator(int n_vertices, int edge_idx) : edge_idx(edge_idx), n_vertices(n_vertices) {}
        edge_iterator &operator++()
        {
            if (edge_idx >= n_vertices * n_vertices)
                return *this;
            ++edge_idx;
            if (edge_idx / n_vertices == edge_idx % n_vertices)
                ++edge_idx;
            return *this;
        }
        int operator*() const {return edge_idx;}
        bool operator==(const edge_iterator &iter) const
        {
            return edge_idx == iter.edge_idx;
        }
    };
};


pair<complete_graph::edge_iterator, complete_graph::edge_iterator>
edges(const complete_graph &g)
{
    return make_pair(complete_graph::edge_iterator(g.n_vertices, 1),
                     complete_graph::edge_iterator(g.n_vertices, g.n_vertices*g.n_vertices));
}

unsigned
num_edges(const complete_graph &g)
{
    return (g.n_vertices - 1) * (g.n_vertices - 1);
}

int
source(int edge, const complete_graph &g)
{
    return edge / g.n_vertices;
}

int
target(int edge, const complete_graph &g)
{
    return edge % g.n_vertices;
}

namespace boost {
    template<>
    struct graph_traits<complete_graph> {
        typedef int                             vertex_descriptor;
        typedef int                             edge_descriptor;
        typedef directed_tag                    directed_category;
        typedef disallow_parallel_edge_tag      edge_parallel_category;
        typedef edge_list_graph_tag             traversal_category;
        typedef complete_graph::edge_iterator   edge_iterator;
        typedef unsigned                        edges_size_type;

        static vertex_descriptor null_vertex() {return -1;}
    };
}



//*****************************************************************************
//                        Function declarations
//*****************************************************************************

/*
 * set_params_from_counts()
 *
 * Sets the parameters of 'hot' according to
 * data-counts. count_u[u][i] should be the number of times u has
 * value i in the data. count_uv[u][v] should be the number of times
 * that u and v have values i and j respecitively in the data. The
 * z-parameters of the tree is set to reflect the data, whereas the
 * x-parameters are set using p_x0_z0 and p_x1_z1.
 */
void set_params_from_counts(Hot &hot,
                            const multi_array<unsigned, 2> &count_u,
                            const multi_array<unsigned, 4> &count_uv,
                            real_type p_x0_z0,
                            real_type p_x1_z1);

/*
 * adjust_params()
 *
 * Adjust the parameters of 'hot' to ensure that the root and its
 * children have their special parameters, and that no probability is
 * zero. Any zero probability found is replaced by 'zero_replacement'.
 */
void adjust_params(Hot &hot, real_type zero_replacement);

/*
 * compute_probs<>()
 *
 * Function used by find_tree that computes the probabilities
 * Pr[Z(u)|X] and Pr[Z(u),Z(v)|X]. The input is the two dimensional
 * matrix containing the data 'D', the dataset to be considered 'd',
 * and the HOT for which the probabilites are to be computed.
 */
template<typename Real, class Matrix>
Real
compute_probs(const Matrix                      &D,
              unsigned                           d,
              const Rooted_tree< Param<Real> >  &hot,
              boost::multi_array<Real, 2>       &p_u_X,
              boost::multi_array<Real, 4>       &p_uv_X);

/*
 * find_tree<>()
 *
 * EM-algorithm for finding a single best tree given weights on the
 * datasets. 'gamma' is a vector/array that holds real numbers with
 * which each dataset is weighted. This is used in the mixture
 * variation for finding mixtures of trees.
 */
template<typename Real>
Real
find_tree(unsigned tree_idx);

/*
 * extend_parent()
 *
 * Assuming that the parent vector represents a tree for the variables
 * that are set in the column mask, the function extends the parent
 * vector and places all other vertices (i.e., those that are zero in
 * column_mask) directly below the root. The parent can then be passed
 * to a rooted tree to get a tree with g_in.num_variables vertices.
 */
void extend_parent(vector<unsigned> &parent,
                   dynamic_bitset<> column_mask,
                   unsigned root);

//*****************************************************************************
//         Template and inline function and member definitions.
//*****************************************************************************

template<typename Real, class Matrix>
Real
compute_probs(const Matrix                      &D,
              unsigned                           d,
              const Rooted_tree< Param<Real> >  &hot,
              boost::multi_array<Real, 2>       &p_u_X,
              boost::multi_array<Real, 4>       &p_uv_X)
{
#ifdef VERBOSE_COMPUTE_PROBS
    clog << "The tree inside compute_probs:\n";
    for (unsigned k = 0; k < hot.size(); ++k)
        clog << "Pr[ Z(" << k << ") | Z(p) ]:\t"
             << hot.data(k).p_z_zp[0][0] << "\t"
             << hot.data(k).p_z_zp[0][1] << "\t"
             << hot.data(k).p_z_zp[1][0] << "\t"
             << hot.data(k).p_z_zp[1][1] << "\n"
             << "Pr[ X(" << k << ") | Z(" << k << ") ]:\t"
             << hot.data(k).p_x_z[0][0] << "\t"
             << hot.data(k).p_x_z[0][1] << "\t"
             << hot.data(k).p_x_z[1][0] << "\t"
             << hot.data(k).p_x_z[1][1] << "\n\n";

    clog << "dataset: ";
    for (unsigned i = 0; i < hot.size(); ++i)
        {
            clog << D[d][i];
        }
    clog << "\n\n";
#endif
    using boost::multi_array;

    const unsigned N = hot.size();
    const unsigned root = hot.root();
    /*
     *  p_Xu_u[u][a]       = Pr[ X_u | Z(u) = a ]
     *   p_X_u[u][a]       = Pr[ X | Z(u) = a ]
     *     p_u[u][a]       = Pr[ Z(u) = a ]
     *   p_v_u[v][u][a][b] = Pr[ Z(v) = a | Z(u) = b ] {u proper ancestor of v}
     * p_Xu_uv[u][v][a][b] = Pr[ X_u | Z(u) = a, Z(v) = b ] {u prop anc of v}
     */
    static multi_array<Real, 2> p_Xu_u(boost::extents[N][2]);
    static multi_array<Real, 2> p_u(boost::extents[N][2]);
    static multi_array<Real, 4> p_v_u(boost::extents[N][N][2][2]);
    static multi_array<Real, 4> p_Xu_uv(boost::extents[N][N][2][2]);
    static multi_array<Real, 2> p_X_u(boost::extents[N][2]);

    /*
     * Compute p_Xu_u, and p_X.
     */
    for (unsigned u = hot.postorder_begin(); u != hot.NONE; u = hot.postorder_next(u))
        {
            for (unsigned a = 0; a < 2; ++a)
                {
                    p_Xu_u[u][a] = hot.data(u).p_x_z[D[d][u]][a];

                    BOOST_FOREACH (unsigned v, hot.children(u))
                        {
                            p_Xu_u[u][a] *=
                                p_Xu_u[v][0] * hot.data(v).p_z_zp[0][a] +
                                p_Xu_u[v][1] * hot.data(v).p_z_zp[1][a];
                        }
                }
        }
    Real p_X =
        hot.data(root).p_z_zp[0][1] * p_Xu_u[root][0] +
        hot.data(root).p_z_zp[1][1] * p_Xu_u[root][1];

#ifdef VERBOSE_COMPUTE_PROBS
    for (unsigned k = 0; k < N; ++k)
        {
            clog << "p_Xu_u[" << k << "]: "
                 << p_Xu_u[k][0] << "\t"
                 << p_Xu_u[k][1] << "\n";
        }
    clog << "Pr[X]: " << p_X << "\n\n";
#endif

    /*
     * Compute p_u.
     */
    p_u[root][0] = hot.data(root).p_z_zp[0][1];
    p_u[root][1] = hot.data(root).p_z_zp[1][1];
    for (unsigned u = hot.preorder_begin(); u != hot.NONE; u = hot.preorder_next(u))
        {
            if (u == root)
                continue;

            unsigned p = hot.parent(u);
            for (unsigned a = 0; a < 2; ++a)
                {
                    p_u[u][a] =
                        hot.data(u).p_z_zp[a][0] * p_u[p][0] +
                        hot.data(u).p_z_zp[a][1] * p_u[p][1];
                }
        }

#ifdef VERBOSE_COMPUTE_PROBS
    for (unsigned k = 0; k < N; ++k)
        {
            clog << "p_u[" << k << "]: "
                 << p_u[k][0] << "\t"
                 << p_u[k][1] << "\n";
        }
    clog << "\n";
#endif

    /*
     * Compute p_v_u {v < u}
     */
    for (unsigned v = 0; v < N; ++v)
        {
            if (v == root)
                continue;

            unsigned u = hot.parent(v);
            p_v_u[v][u][0][0] = hot.data(v).p_z_zp[0][0];
            p_v_u[v][u][0][1] = hot.data(v).p_z_zp[0][1];
            p_v_u[v][u][1][0] = hot.data(v).p_z_zp[1][0];
            p_v_u[v][u][1][1] = hot.data(v).p_z_zp[1][1];

            unsigned w = u;
            for (u = hot.parent(u); u != hot.NONE; u = hot.parent(u))
                {
                    const Param<Real> &wparam = hot.data(w);
                    for (unsigned a = 0; a < 2; ++a)
                        for (unsigned b = 0; b < 2; ++b)
                            p_v_u[v][u][a][b] =
                                p_v_u[v][w][a][0] * wparam.p_z_zp[0][b] +
                                p_v_u[v][w][a][1] * wparam.p_z_zp[1][b];
                    w = u;
                }
        }

#ifdef VERBOSE_COMPUTE_PROBS
    for (unsigned v = 0; v < N; ++v)
        for (unsigned u = 0; u < N; ++u)
            {
                clog << "p_v_u[" << v << "][" << u << "]:    ";
                for (unsigned a = 0; a < 2; ++a)
                    for (unsigned b = 0; b < 2; ++b)
                    {
                        clog << p_v_u[v][u][a][b] << "\t";
                    }
                clog << "\n";
            }
    clog << "\n";
#endif

    /*
     * Compute p_Xu_uv {v < u}
     */
    for (unsigned v = hot.postorder_begin(); v != hot.NONE; v = hot.postorder_next(v))
        for (unsigned u = hot.parent(v); u != hot.NONE; u = hot.parent(u))
            {
                unsigned xu = D[d][u];
                for (unsigned a = 0; a < 2; ++a)
                    for (unsigned b = 0; b < 2; ++b)
                        {
                            if (u == root && a == 0)
                                {
                                    p_Xu_uv[u][v][a][b] = Real(0);
                                    continue;
                                }

                            p_Xu_uv[u][v][a][b] = hot.data(u).p_x_z[xu][a];
                            BOOST_FOREACH (unsigned child, hot.children(u))
                                if (child == v)
                                    {
                                        p_Xu_uv[u][v][a][b] *=
                                            p_Xu_u[v][b];
                                    }
                                else if (hot.descendant(v, child))
                                    {
                                        p_Xu_uv[u][v][a][b] *=
                                            p_Xu_uv[child][v][1][b] *
                                            p_v_u[v][child][b][1] *
                                            p_v_u[child][u][1][a] /
                                            p_v_u[v][u][b][a]
                                            +
                                            p_Xu_uv[child][v][0][b] *
                                            p_v_u[v][child][b][0] *
                                            p_v_u[child][u][0][a] /
                                            p_v_u[v][u][b][a];
                                    }
                                else
                                    {
                                        p_Xu_uv[u][v][a][b] *=
                                            p_Xu_u[child][0] *
                                            p_v_u[child][u][0][a]
                                            +
                                            p_Xu_u[child][1] *
                                            p_v_u[child][u][1][a];
                                    }
                        }
            }

#ifdef VERBOSE_COMPUTE_PROBS
    for (unsigned u = 0; u < N; ++u)
        for (unsigned v = 0; v < N; ++v)
            {
                clog << "p_Xu_uv[" << u << "][" << v << "]:  ";
                for (unsigned a = 0; a < 2; ++a)
                    for (unsigned b = 0; b < 2; ++b)
                    {
                        clog << p_Xu_uv[u][v][a][b] << "\t";
                    }
                clog << "\n";
            }
    clog << "\n";
#endif

    /*
     * Compute p_X_u.
     */
    for (unsigned u = 0; u < N; ++u)
        {
            if (u == root)
                continue;
            for (unsigned a = 0; a < 2; ++a)
                {
                    p_X_u[u][a] =
                        (p_Xu_uv[root][u][0][a] * p_v_u[u][root][a][0]
                         * p_u[root][0] +
                         p_Xu_uv[root][u][1][a] * p_v_u[u][root][a][1]
                         * p_u[root][1]);
                    p_X_u[u][a] /= p_u[u][a];
                }
        }
    p_X_u[root][0] = p_Xu_u[root][0];
    p_X_u[root][1] = p_Xu_u[root][1];

    /*
     * Compute return values p_u_X.
     */
    for (unsigned u = 0; u < N; ++u)
        {
            if (u == root)
                continue;

            for (unsigned a = 0; a < 2; ++a)
                {
                    p_u_X[u][a] = p_X_u[u][a] * p_u[u][a] / p_X;
                }
        }
    p_u_X[root][0] = p_X_u[root][0] * hot.data(root).p_z_zp[0][1] / p_X;
    p_u_X[root][1] = p_X_u[root][1] * hot.data(root).p_z_zp[1][1] / p_X;

#ifdef VERBOSE_COMPUTE_PROBS
    for (unsigned k = 0; k < N; ++k)
        {
            std::clog << "p_u_X[" << k << "]: "
                      << p_u_X[k][0] << "\t"
                      << p_u_X[k][1] << "\n";
        }
    clog << "\n";
#endif

    /*
     * Now to the grizzly part. compute return values p_uv_X.
     */
    for (unsigned u = 0; u < N; ++u)
        {
            p_uv_X[u][u][0][0] = Real(0.0);
            p_uv_X[u][u][0][1] = Real(0.0);
            p_uv_X[u][u][1][0] = Real(0.0);
            p_uv_X[u][u][1][1] = Real(0.0);
        }

    for (unsigned v = 0; v < N; ++v)
        {
            for (unsigned u = hot.parent(v); u != hot.NONE; u = hot.parent(u))
                {
                    if (u == root)
                        {
                            for (unsigned a = 0; a < 2; ++a)
                                for (unsigned b = 0; b < 2; ++b)
                                    {
                                        p_uv_X[u][v][a][b] =
                                            p_Xu_uv[u][v][a][b] *
                                            p_v_u[v][u][b][a] *
                                            p_u[u][a] /
                                            p_X;
                                    }

//                             p_uv_X[u][v][0][0] = Real(0.0);
//                             p_uv_X[u][v][0][1] = Real(0.0);
//                             p_uv_X[u][v][1][0] =
//                                 p_Xu_uv[u][v][1][0] *
//                                 p_v_u[v][u][0][1] *
//                                 p_u[u][1] /
//                                 p_X;
//                             p_uv_X[u][v][1][1] =
//                                 p_Xu_uv[u][v][1][1] *
//                                 p_v_u[v][u][1][1] *
//                                 p_u[u][1] /
//                                 p_X;
                        }
                    else
                        {
                            for (unsigned a = 0; a < 2; ++a)
                                for (unsigned b = 0; b < 2; ++b)
                                    {
                                        Real p_X_uv =
                                            p_X_u[u][a] * p_Xu_uv[u][v][a][b] /
                                            p_Xu_u[u][a];
                                        Real uv =
                                            p_v_u[v][u][b][a] *
                                            p_u[u][a];
                                        p_uv_X[u][v][a][b] =
                                            p_X_uv * uv / p_X;
                                    }
                         }

                    p_uv_X[v][u][0][0] = p_uv_X[u][v][0][0];
                    p_uv_X[v][u][0][1] = p_uv_X[u][v][1][0];
                    p_uv_X[v][u][1][0] = p_uv_X[u][v][0][1];
                    p_uv_X[v][u][1][1] = p_uv_X[u][v][1][1];
                }
        }

    for (unsigned u = 0; u < N-1; ++u)
        for (unsigned v = u+1; v < N; ++v)
            {
                unsigned w = hot.lca(u, v);
                if (u == w || v == w)
                    continue;

                for (unsigned a = 0; a < 2; ++a)
                    for (unsigned b = 0; b < 2; ++b)
                        {
                            p_uv_X[u][v][a][b] =
                                p_uv_X[w][u][1][a] *
                                p_uv_X[w][v][1][b] /
                                p_u_X[w][1];
                            if (w != root)
                                {
                                    p_uv_X[u][v][a][b] +=
                                        p_uv_X[w][u][0][a] *
                                        p_uv_X[w][v][0][b] /
                                        p_u_X[w][0];
                                }
                            p_uv_X[v][u][b][a] = p_uv_X[u][v][a][b];
                        }
            }

#ifdef VERBOSE_COMPUTE_PROBS
    for (unsigned u = 0; u < N; ++u)
        for (unsigned v = 0; v < N; ++v)
            {
                clog << "p_uv_X[" << u << "][" << v << "]:    ";
                for (unsigned a = 0; a < 2; ++a)
                    for (unsigned b = 0; b < 2; ++b)
                    {
                        clog << p_uv_X[u][v][a][b] << "\t";
                    }
                clog << "\n";
            }
    clog << "\n";
#endif
    return p_X;
}



template<typename Real>
Real
find_tree(unsigned tree_idx)
{
    using boost::multi_array;

    /* A few variables we will be using throughout the function. */
    const unsigned       N = g_in.num_variables;
    const unsigned       zero_one[2] = {0, 1};
    Hot                 &hot = g_trees[tree_idx];
    /*
     *   p_u_X[u][a]        =  Pr[ Z(u)=a         | X ]
     *  p_uv_X[u][v][a][b]  =  Pr[ Z(u)=a, Z(v)=b | X ]
     * theta_z[u][v][a][b]  =  Pr[ Z(v)=b         | Z(u)=a ]
     * theta_x[u][s][a]     =  Pr[ X(u)=s         | Z(u)=a ]
     *       A[u][v][a][b]  =  SUM_i gamma[i] * Pr[ Z(u)=a, Z(v)=b | X_i ]
     *       B[u][s][a]     =  SUM_i gamma[i] * Pr[ Z(u)=a | X_i ] {X_i(u)=s}
     */
    multi_array<Real, 2> p_u_X(boost::extents[N][2]);
    multi_array<Real, 4> p_uv_X(boost::extents[N][N][2][2]);
    multi_array<Real, 4> A(boost::extents[N][N][2][2]);
    multi_array<Real, 3> B(boost::extents[N][2][2]);
    multi_array<Real, 4> theta_z(boost::extents[N][N][2][2]);
    multi_array<Real, 3> theta_x(boost::extents[N][2][2]);
    multi_array<Real, 2> weights(boost::extents[N][N]);
    std::vector<unsigned> parent(N);

    /* p_X_T[x][i] == Pr[ data[x] | T_i ]*/
    multi_array<Real, 2> p_X_T(boost::extents[g_in.num_datasets][g_in.num_trees]);

    /* Compute p_X_T for the all trees. */
    for (unsigned x = 0; x < g_in.num_datasets; ++x)
    {
        for (unsigned i = 0; i < g_in.num_trees; ++i)
        {
            p_X_T[x][i] = compute_p_X(g_trees[i], g_in.data, x);
        }
    }

    Real new_L(1.0);
    Real old_L(0.0);

    for (unsigned iterations = 0; ; ++iterations)
        {
            std::clog << "starting iteration " << iterations << "\n\n";

            clog << "parent structure of the tree:\n";
            for (unsigned k = 0; k < N; ++k)
                if (k != g_root)
                    std::clog << setw(3) << k;
            std::clog << "\n";
            for (unsigned k = 0; k < N; ++k)
                if (k != g_root)
                    std::clog << setw(3) << hot.parent(k);
            std::clog << "\n";

//             for (unsigned k = 0; k < N; ++k)
//                 clog << "Pr[ Z(" << k << ") | Z(p) ]:\t"
//                      << hot.data(k).p_z_zp[0][0] << "\t"
//                      << hot.data(k).p_z_zp[0][1] << "\t"
//                      << hot.data(k).p_z_zp[1][0] << "\t"
//                      << hot.data(k).p_z_zp[1][1] << "\n";
//             clog << "\n";
//             for (unsigned k = 0; k < N; ++k)
//                 clog << "Pr[ X(" << k << ") | Z(" << k << ") ]:\t"
//                      << hot.data(k).p_x_z[0][0] << "\t"
//                      << hot.data(k).p_x_z[0][1] << "\t"
//                      << hot.data(k).p_x_z[1][0] << "\t"
//                      << hot.data(k).p_x_z[1][1] << "\n";
//             clog << "\n";


            /* Reset A and B. */
            std::fill(A.origin(), A.origin() + N*N*4, Real(0.0));
            std::fill(B.origin(), B.origin() + N*4, Real(0.0));

            /*
             * Compute Pr[Z(u) | X] and Pr[Z(u), Z(v) | X] for each X,
             * and compute the sums in A and B. Also, set new_L to
             * the likelihood of the HOT, i.e, new_L <- Pr[D].
             */
            Real p_X(0.0);
            for (unsigned i = 0; i < g_in.num_datasets; ++i)
                {
                    if (g_gamma[tree_idx][i] == Real(0))
                        continue;
                    p_X = compute_probs(g_in.data, i, hot, p_u_X, p_uv_X);
                    p_X_T[i][tree_idx] = p_X;

#ifdef VERBOSE_FIND_TREE
                    clog << "data: ";
                    for (unsigned u = 0; u < N; ++u)
                        {
                            clog << g_in.data[i][u];
                        }
                    clog << "   " << g_in.input_freq[i] << "\n\n";

                    for (unsigned u = 0; u < N; ++u)
                        {
                            clog << "Pr[Z(" << u << ")=a|X]:\t"
                                 << p_u_X[u][0] << "\t" << p_u_X[u][1] << "\n";
                        }
                    clog << "\n";

                    for (unsigned u = 0; u < N-1; ++u)
                        {
                            for (unsigned v = u+1; v < N; ++v)
                                {
                                    clog << "Pr[Z(" << u << ")=a, Z("
                                         << v << ")=b|X]:\t"
                                         << p_uv_X[u][v][0][0] << "\t"
                                         << p_uv_X[u][v][0][1] << "\t"
                                         << p_uv_X[u][v][1][0] << "\t"
                                         << p_uv_X[u][v][1][1] << "\n";
                                }
                        }
                    clog << "\n";
#endif

                    for (unsigned u = 0; u < N-1; ++u)
                        for (unsigned v = u+1; v < N; ++v)
                            {
                                BOOST_FOREACH (unsigned a, zero_one)
                                    {
                                        BOOST_FOREACH (unsigned b, zero_one)
                                            {
                                                A[u][v][a][b] +=
                                                    g_in.input_freq[i]*g_gamma[tree_idx][i]*p_uv_X[u][v][a][b];
                                                A[v][u][b][a] = A[u][v][a][b];
                                            }
                                    }
                            }

                    for (unsigned u = 0; u < N; ++u)
                        BOOST_FOREACH (unsigned a, zero_one)
                            {
                                B[u][ g_in.data[i][u] ][a] +=
                                    g_in.input_freq[i]*g_gamma[tree_idx][i]*p_u_X[u][a];
                            }
                }

#ifdef VERBOSE_FIND_TREE
            for (unsigned u = 0; u < N; ++u)
                for (unsigned v = 0; v < N; ++v)
                    {
                        clog << "A[" << u << "][" << v << "]: "
                             << A[u][v][0][0] << "\t"
                             << A[u][v][0][1] << "\t"
                             << A[u][v][1][0] << "\t"
                             << A[u][v][1][1] << "\n";
                    }
            clog << "\n";

            for (unsigned u = 0; u < N; ++u)
                {
                    clog << "B[" << u << "]:   "
                         << B[u][0][0] << "\t"
                         << B[u][0][1] << "\t"
                         << B[u][1][0] << "\t"
                         << B[u][1][1] << "\n";
                }
            clog << "\n";
#endif

            new_L = Real(1.0);
            for (unsigned x = 0; x < g_in.num_datasets; ++x)
            {
                Real p(0.0);
                for (unsigned i = 0; i < g_in.num_trees; ++i)
                {
                    p += g_lambda[i] * p_X_T[x][i];
                }
                for (unsigned i = 0; i < g_in.input_freq[x]; ++i)
                {
                    new_L *= p;
                }
            }

            clog << "========== old likelihood: " << old_L << "\t"
                 << "new likelihood: " << new_L << " ==========\n\n";

            /* Break the loop if cutoff reached. */
            if (iterations >= g_in.max_tree_iterations)
                break;
            Real diff;
            if (new_L < old_L)
                {
                    clog << "========== warning: likelihood decreased! ==========\n";
                    diff = old_L - new_L;
                }
            else
                {
                    diff = new_L - old_L;
                }
            if (iterations > 0 && diff/old_L < Real(g_in.tree_cutoff))
                break;

            /* Save old likelihood. */
            old_L = new_L;

            /* Get new weights for all possible arcs and save the
               quantities that correspond to parameters on arcs for
               later use. */

            Real pseudo_count = Real(g_in.pseudo_count);
            for (unsigned u = 0; u < N; ++u)
                {
                    for (unsigned v = 0; v < N; ++v)
                        {
                            A[u][v][0][0] += pseudo_count;
                            A[u][v][0][1] += pseudo_count;
                            A[u][v][1][0] += pseudo_count;
                            A[u][v][1][1] += pseudo_count;
                        }
                }
            for (unsigned u = 0; u < N; ++u)
                {
                    B[u][0][0] += pseudo_count;
                    B[u][0][1] += pseudo_count;
                    B[u][1][0] += pseudo_count;
                    B[u][1][1] += pseudo_count;
                }

            for (unsigned u = 0; u < N; ++u)
                for (unsigned v = 0; v < N; ++v)
                    {
                        if (u == v)
                            continue;

                        if (u == g_root)
                            {
                                Real sum = B[v][0][0] + B[v][1][0] +
                                    B[v][0][1] + B[v][1][1];
                                theta_z[u][v][1][0] =
                                    (B[v][0][0] + B[v][1][0]) / sum;
                                theta_z[u][v][1][1] =
                                    Real(1.0) - theta_z[u][v][1][0];
                                theta_z[u][v][0][0] = Real(0);
                                theta_z[u][v][0][1] = Real(0);
                                continue;
                            }

                        Real sum0 = A[u][v][0][0] + A[u][v][0][1];
                        Real sum1 = A[u][v][1][0] + A[u][v][1][1];

                        if (sum0 == Real(0))
                        {
                            sum0 = sum1 / Real(g_in.num_datasets);
                        }
                        if (sum1 == Real(0))
                        {
                            sum1 = sum0 / Real(g_in.num_datasets);
                        }

                        theta_z[u][v][0][1] = A[u][v][0][1] / sum0;
                        theta_z[u][v][0][0] = A[u][v][0][0] / sum0;
                        theta_z[u][v][1][0] = A[u][v][1][0] / sum1;
                        theta_z[u][v][1][1] = A[u][v][1][1] / sum1;

                        if (Real(g_in.p_z1_zp0) >= Real(0) &&
                            theta_z[u][v][0][1] > Real(g_in.p_z1_zp0))
                            {
                                theta_z[u][v][0][1] = g_in.p_z1_zp0;
                                theta_z[u][v][0][0] = Real(1.0) - Real(g_in.p_z1_zp0);
                            }
                    }

            theta_x[g_root][0][0] = Real(1);
            theta_x[g_root][0][1] = Real(0);
            theta_x[g_root][1][0] = Real(0);
            theta_x[g_root][1][1] = Real(1);

            if (g_in.global_x_params)
                {
                    Real sum00(0.0);
                    Real sum01(0.0);
                    Real sum10(0.0);
                    Real sum11(0.0);

                    for (unsigned u = 0; u < N; ++u)
                        {
                            if (u == g_root)
                                continue;

                            sum00 += B[u][0][0];
                            sum01 += B[u][0][1];
                            sum10 += B[u][1][0];
                            sum11 += B[u][1][1];
                        }

                    for (unsigned u = 0; u < N; ++u)
                        {
                            if (u == g_root)
                                continue;

                            theta_x[u][0][0] = sum00 / (sum00 + sum10);
                            theta_x[u][1][0] = sum10 / (sum00 + sum10);
                            theta_x[u][0][1] = sum01 / (sum01 + sum11);
                            theta_x[u][1][1] = sum11 / (sum01 + sum11);

                            if (Real(g_in.p_x1_z0) >= Real(0) &&
                                theta_x[u][1][0] > Real(g_in.p_x1_z0))
                                {
                                    theta_x[u][1][0] = g_in.p_x1_z0;
                                    theta_x[u][0][0] = Real(1.0) - Real(g_in.p_x1_z0);
                                }
                        }
                }
            else
                {
                    for (unsigned u = 0; u < N; ++u)
                        {
                            if (u == g_root)
                                continue;

                            Real sum0 = B[u][0][0] + B[u][1][0];
                            Real sum1 = B[u][0][1] + B[u][1][1];
                            theta_x[u][0][0] = B[u][0][0] / sum0;
                            theta_x[u][1][0] = B[u][1][0] / sum0;
                            theta_x[u][0][1] = B[u][0][1] / sum1;
                            theta_x[u][1][1] = B[u][1][1] / sum1;

                            if (Real(g_in.p_x1_z0) >= Real(0)
                                && theta_x[u][1][0] > Real(g_in.p_x1_z0))
                                {
                                    theta_x[u][1][0] = Real(g_in.p_x1_z0);
                                    theta_x[u][0][0] = Real(1.0) - Real(g_in.p_x1_z0);
                                }
                        }
                }

#ifdef VERBOSE_FIND_TREE
            for (unsigned u = 0; u < N; ++u)
                for (unsigned v = 0; v < N; ++v)
                    {
                        clog << "theta_z[" << u << "][" << v << "] : "
                             << theta_z[u][v][0][0] << " "
                             << theta_z[u][v][0][1] << " "
                             << theta_z[u][v][1][0] << " "
                             << theta_z[u][v][1][1] << "\n";
                    }
            clog << "\n";

            for (unsigned u = 0; u < N; ++u)
                {
                    clog << "theta_x[" << u << "] : "
                         << theta_x[u][0][0] << " "
                         << theta_x[u][0][1] << " "
                         << theta_x[u][1][0] << " "
                         << theta_x[u][1][1] << "\n";
                }
            clog << "\n";
#endif

            for (unsigned u = 0; u < N; ++u)
                {
                    for (unsigned v = 0; v < N; ++v)
                        {
                            if (u == v)
                                continue;

                            if (v == g_root)
                                {
                                    weights[u][v] = Real(0.0);
                                    continue;
                                }

                            if (u == g_root)
                                {
                                    weights[u][v] =
                                        (B[v][0][0] + B[v][1][0]) * log(theta_z[u][v][1][0]) +
                                        (B[v][0][1] + B[v][1][1]) * log(theta_z[u][v][1][1]);
//                                     weights[u][v] +=
//                                         B[v][0][0] * log(theta_x[v][0][0]) +
//                                         B[v][0][1] * log(theta_x[v][0][1]) +
//                                         B[v][1][0] * log(theta_x[v][1][0]) +
//                                         B[v][1][1] * log(theta_x[v][1][1]);
                                    continue;
                                }

                            weights[u][v] =
                                A[u][v][0][0] * log(theta_z[u][v][0][0]) +
                                A[u][v][0][1] * log(theta_z[u][v][0][1]) +
                                A[u][v][1][0] * log(theta_z[u][v][1][0]) +
                                A[u][v][1][1] * log(theta_z[u][v][1][1]);
//                             weights[u][v] +=
//                                 B[v][0][0] * log(theta_x[v][0][0]) +
//                                 B[v][0][1] * log(theta_x[v][0][1]) +
//                                 B[v][1][0] * log(theta_x[v][1][0]) +
//                                 B[v][1][1] * log(theta_x[v][1][1]);
                        }
                }

#ifdef VERBOSE_FIND_TREE
            for (unsigned u = 0; u < N; ++u)
                {
                    for (unsigned v = 0; v < N; ++v)
                        {
                            std::clog << weights[u][v] << "\t";
                        }
                    clog << "\n";
                }
            clog << "\n";
#endif

            /* Get a new rooted tree based on the weights and set the
               parameters of the arcs appropriately (from saved
               value). */
            if (!g_in.const_tree)
                {
                    typedef boost::graph_traits<complete_graph>::edge_descriptor Edge;
                    typedef boost::graph_traits<complete_graph>::vertex_descriptor Vertex;

                    const dynamic_bitset<> &mask = g_in.column_masks[tree_idx];
                    if ((~mask).none())
                        {
                            complete_graph g(N);
                            Vertex branching_roots[] = {g_root};
                            vector<Edge> branching_edges;
                            edmonds_optimum_branching<true, true, true>
                                (g, boost::identity_property_map(), weights.origin(),
                                 branching_roots, branching_roots + 1, back_inserter(branching_edges));
                            BOOST_FOREACH (Edge e, branching_edges)
                            {
                                parent[target(e, g)] = source(e, g);
                            }
                            //                            max_spanning_arborescence<Real>(weights, N, g_root, parent);
                        }
                    else
                        {
                            /* get new weights corresponding to the column mask. */
                            unsigned small_size = mask.count();
                            multi_array<Real, 2> new_weights(boost::extents[small_size][small_size]);
                            vector<unsigned> indices;
                            for (dynamic_bitset<>::size_type i = mask.find_first(); i != dynamic_bitset<>::npos;
                                 i = mask.find_next(i))
                                {
                                    indices.push_back(i);
                                }
                            unsigned row = 0;
                            BOOST_FOREACH (unsigned row_idx, indices)
                                {
                                    unsigned col = 0;
                                    BOOST_FOREACH (unsigned col_idx, indices)
                                        {
                                            new_weights[row][col] = weights[row_idx][col_idx];
                                            col++;
                                        }
                                    row++;
                                }

                            complete_graph g(small_size);
                            Vertex branching_roots[] = {g_root};
                            vector<Edge> branching_edges;
                            edmonds_optimum_branching<true, true, true>
                                (g, boost::identity_property_map(), new_weights.origin(),
                                 branching_roots, branching_roots + 1, back_inserter(branching_edges));
                            BOOST_FOREACH (Edge e, branching_edges)
                            {
                                parent[target(e, g)] = source(e, g);
                            }
                            //                            max_spanning_arborescence<Real>(new_weights, small_size, g_root, parent);
                            extend_parent(parent, mask, g_root);
                        }
                    hot.reset(parent, g_root);
                }
            for (unsigned u = 0; u < N; ++u)
                {
                    Param<Real> &param = hot.data(u);

                    if (u == g_root)
                        {
                            param.p_z_zp[0][0] = Real(0.0);
                            param.p_z_zp[1][0] = Real(0.0);
                            param.p_z_zp[0][1] = Real(0);
                            param.p_z_zp[1][1] = Real(1);

                            param.p_x_z[0][0] = Real(1);
                            param.p_x_z[0][1] = Real(0);
                            param.p_x_z[1][0] = Real(0);
                            param.p_x_z[1][1] = Real(1);
                            continue;
                        }

                    unsigned p = hot.parent(u);

                    if (p == g_root)
                        {
                            param.p_z_zp[0][0] = Real(0);
                            param.p_z_zp[0][1] = theta_z[p][u][1][0];
                            param.p_z_zp[1][0] = Real(0);
                            param.p_z_zp[1][1] = theta_z[p][u][1][1];
                            param.p_x_z[0][0] = theta_x[u][0][0];
                            param.p_x_z[0][1] = theta_x[u][0][1];
                            param.p_x_z[1][0] = theta_x[u][1][0];
                            param.p_x_z[1][1] = theta_x[u][1][1];
                            continue;
                        }

                    param.p_z_zp[0][0] = theta_z[p][u][0][0];
                    param.p_z_zp[0][1] = theta_z[p][u][1][0];
                    param.p_z_zp[1][0] = theta_z[p][u][0][1];
                    param.p_z_zp[1][1] = theta_z[p][u][1][1];
                    param.p_x_z[0][0] = theta_x[u][0][0];
                    param.p_x_z[0][1] = theta_x[u][0][1];
                    param.p_x_z[1][0] = theta_x[u][1][0];
                    param.p_x_z[1][1] = theta_x[u][1][1];
                }
            adjust_params(hot, real_type(1) / real_type(g_in.total_datasets * 10));
        }

    return new_L;
}

//*****************************************************************************
//                                main()
//*****************************************************************************

int
main(int argc, char *argv[])
{
    unsigned seed;
    /*
     * Parse the command line options.
     */

    po::options_description visible_opts("Command Line Options");
    po::options_description hidden_opts("");
    po::options_description all_options("");

    try
        {
            /* Declare visible options */
            visible_opts.add_options()
                ("help,h", "display this help and exit")
                ("max-mix-iters",
                 po::value<int>(&g_in.max_mixture_iterations)
                 ->default_value(1000),
                 "Maximum number of mixture iterations")
                ("max-tree-iters",
                 po::value<int>(&g_in.max_tree_iterations)
                 ->default_value(1000),
                 "Maximum number of tree iterations")
                ("mixture-cutoff",
                 po::value<double>(&g_in.mixture_cutoff)
                 ->default_value(1e-3),
                 "Cutoff value for the mixtures")
                ("tree-cutoff",
                 po::value<double>(&g_in.tree_cutoff)
                 ->default_value(1e-2),
                 "Cutoff value for the trees")
                ("global-x-params",
                 po::bool_switch(&g_in.global_x_params),
                 "if set, Pr[X(u)|Z(u)=0] and Pr[X(u)|Z(u)=1] are global, i.e., every edge uses the same parameters.")
                ("pseudo-count",
                 po::value<double>(&g_in.pseudo_count)->default_value(0),
                 "Pseudocounts added to A and B")
                ("p-z1-zp0",
                 po::value<double>(&g_in.p_z1_zp0)->default_value(0.5),
                 "Maximum value for Pr[Z(u)=1|Z(p(u))=0]. Set to a negative number to estimate it from data.")
                ("p-x1-z0",
                 po::value<double>(&g_in.p_x1_z0)->default_value(0.5),
                 "Maximum value for Pr[X(u)=1|Z(u)=0]. Set to a negative number to estimate it from data.")
                ("start-tree", po::value< vector<string> >(&g_in.start_trees),
                 "the type of tree to start the EM-algorithm with. Possible values are 'random', and 'starz'. Any other string is taken to be a filename containing a start tree. This option can be given as many times as NUM_TREES.")
                ("lambda", po::value< vector<double> >(&g_in.start_lambdas),
                 "You can give the prior probability of the trees by giving lambda several times.")
                ("column-mask,c", po::value< vector< dynamic_bitset<> > >(&g_in.column_masks),
                 "The argument is a string of 1:s and 0:s, the length of which must match the number of columns in the data. "
                 "Only the columns with a 1 will be used to build the tree. The vertices representing columns with a zero will all be children of the root. "
                 "When using mixtures of trees, you can give this option as many times as the number of trees.")
                ("const-tree",
                 po::bool_switch(&g_in.const_tree),
                 "If set, the program will not change the start tree but only maximize the parameters of that tree.")
                ("hard-em",
                 po::bool_switch(),
                 "If set, the program will use hard EM to cluster the input instead of the usual soft EM.")
                ("seed", po::value<unsigned>(&seed),
                 "seed for the random number generator")
                ;

            /* Declare positional options. */
            hidden_opts.add_options()
                ("input-filename",po::value<string>())
                ("num-trees", po::value<unsigned>(&g_in.num_trees))
                ;

            po::positional_options_description positional_options;
            positional_options.add("input-filename", 1);
            positional_options.add("num-trees", 1);

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
    if (g_options.count("input-filename") == 0 ||
        g_options.count("num-trees") == 0)
        {
            cerr << PROG_NAME << ": Too few arguments.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Check that the small probabilities were given reasonable values*/
    if (g_in.p_z1_zp0 > 1)
        {
            cerr << PROG_NAME << ": bad value given for --p-z1-zp0.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }
    if (g_in.p_x1_z0 > 1)
        {
            cerr << PROG_NAME << ": bad value given for --p-x1-z0.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* Check num_trees */
    if (g_in.num_trees <= 0)
        {
            cerr << PROG_NAME << ": NUM_TREES must be at least one.\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* check the number of start trees */
    if (g_in.start_trees.size() > g_in.num_trees)
        {
            cerr << PROG_NAME << ": too many start trees given\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    while (g_in.start_trees.size() < g_in.num_trees)
        {
            g_in.start_trees.push_back("random");
        }

    if (g_in.column_masks.size() > g_in.num_trees)
        {
            cerr << PROG_NAME << ": too many column masks given\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    if (g_in.start_lambdas.size() > g_in.num_trees)
        {
            cerr << PROG_NAME << ": too many lambdas given\n"
                 << "Try '" << PROG_NAME << " --help' for more information.\n";
            exit(EXIT_FAILURE);
        }

    /* check that start lambdas are reasonable. */

    const double tolerance = 1e-6;
    double lambda_sum = accumulate(g_in.start_lambdas.begin(),
                                   g_in.start_lambdas.end(),
                                   0.0);

    if (g_in.start_lambdas.size() == g_in.num_trees)
        {
            if (fabs(1.0 - lambda_sum) > tolerance)
                {
                    cerr << PROG_NAME << ": lambdas do not sum to 1.0\n"
                         << "Try '" << PROG_NAME << " --help' for more information.\n";
                    exit(EXIT_FAILURE);
                }

            g_in.start_lambdas.back() =
                1.0 - accumulate(g_in.start_lambdas.begin(),
                                 g_in.start_lambdas.end() - 1,
                                 0.0);
        }
    else
        {
            int missing_lambdas = g_in.num_trees - g_in.start_lambdas.size();
            if (lambda_sum > 1.0 - tolerance * missing_lambdas)
                {
                    cerr << PROG_NAME << ": lambdas too large\n"
                         << "Try '" << PROG_NAME << " --help' for more information.\n";
                    exit(EXIT_FAILURE);
                }

            for (int i = 0; i < missing_lambdas; ++i)
                {
                    g_in.start_lambdas.push_back((1.0 - lambda_sum) / missing_lambdas);
                }
        }

    /* Seed the random number generators. */
    if (g_options.count("seed"))
        {
            init_rand(g_generator, seed);
        }
    else
        {
            init_rand(g_generator);
        }
    /*
     * Read the input file.
     */

    /* Open the file. */
    string filename = g_options["input-filename"].as<string>();
    ifstream infile(filename.c_str());
    if (!infile)
        {
            cerr << PROG_NAME << ": Could not open file '" << filename << "'\n";
            exit(EXIT_FAILURE);
        }

    /* Read all input lines into a vector. */
    vector<string> input_lines;
    read_input_lines(infile, input_lines);
    infile.close();

    /* Do we have at least one input line? */
    if (input_lines.empty())
        {
            cerr << PROG_NAME << ": no valid input lines detected.\n";
            exit(EXIT_FAILURE);
        }

    /* Determine the number of variables and the size of the dataset. */
    g_in.num_variables = input_lines[0].size() + 1;

    /* Do all lines have the same number of variables? */
    BOOST_FOREACH (string s, input_lines)
        {
            if (s.size() != input_lines.front().size())
                {
                    cerr << PROG_NAME << ": error reading input file\n";
                    cerr << PROG_NAME << ": every valid input line must contain the same number of variables\n";
                    exit(EXIT_FAILURE);
                }
        }

    /* sort and uniquify the input */
    sort(input_lines.begin(), input_lines.end());
    vector<string> unique_lines;
    unique_copy(input_lines.begin(), input_lines.end(),
                back_inserter(unique_lines));
    g_in.input_freq.resize(unique_lines.size());
    unsigned j = 0;
    for (unsigned i = 0; i < input_lines.size(); )
        {
            g_in.input_freq[j] += 1;
            ++i;
            while (i < input_lines.size() && input_lines[i] == input_lines[i-1])
                {
                    g_in.input_freq[j] += 1;
                    ++i;
                }
            ++j;
        }
    g_in.num_datasets = unique_lines.size();
    g_in.total_datasets = input_lines.size();

    /* Create the input matrix. */
    g_in.data.resize(boost::extents[unique_lines.size()][g_in.num_variables]);
    unsigned r = 0;
    BOOST_FOREACH (string s, unique_lines)
        {
            unsigned c = 0;
            g_in.data[r][c] = true;
            ++c;
            BOOST_FOREACH (char ch, s)
                {
                    g_in.data[r][c] = ch == '0' ? false : true;
                    ++c;
                }
            ++r;
        }

    clog << "Read " << input_lines.size() << " datasets "
         << "with " << g_in.num_variables << " vertices.\n";

    /* Verify the column masks. */
    BOOST_FOREACH (dynamic_bitset<> &mask, g_in.column_masks)
        {
            if (mask.size() != g_in.num_variables - 1)
                {
                    cerr << PROG_NAME << ": column mask with bad size\n"
                         << "Try '" << PROG_NAME << " --help' for more information.\n";
                    exit(EXIT_FAILURE);
                }
        }
    while (g_in.column_masks.size() < g_in.num_trees)
        {
            g_in.column_masks.push_back(~dynamic_bitset<>(g_in.num_variables - 1));
        }
    BOOST_FOREACH (dynamic_bitset<> &mask, g_in.column_masks)
        {
            mask.push_back(1);
            dynamic_bitset<> mask_copy = mask;
            for (unsigned i = 0; i < mask.size(); ++i)
                {
                    mask[i] = mask_copy[mask.size() - i - 1];
                }
        }


    /* Count the the combinations in the input for each pair of vertices. */
    multi_array<unsigned, 4> count_uv(boost::extents[g_in.num_variables][g_in.num_variables][2][2]);
    multi_array<unsigned, 2> count_u(boost::extents[g_in.num_variables][2]);

    for (unsigned i = 0; i < g_in.num_datasets; ++i)
        {
            for (unsigned u = 0; u < g_in.num_variables; ++u)
                {
                    count_u[u][g_in.data[i][u]] += g_in.input_freq[i];
                    for (unsigned v = u+1; v < g_in.num_variables; ++v)
                        {
                            count_uv[u][v][g_in.data[i][u]][g_in.data[i][v]] += g_in.input_freq[i];
                            count_uv[v][u][g_in.data[i][v]][g_in.data[i][u]] += g_in.input_freq[i];
                        }
                }
        }

    /* Get starting trees. */
    g_trees.resize(g_in.num_trees);
    for (unsigned i = 0; i < g_in.num_trees; ++i)
        {
            if (g_in.start_trees[i] == "random")
                {
                    vector<Hot::vid_t> parent;
                    gen_random_tree<real_type>(g_in.column_masks[i].count(),
                                               g_root,
                                               parent);
                    extend_parent(parent, g_in.column_masks[i], g_root);
                    g_trees[i].reset(parent, g_root);
                    set_params_from_counts(g_trees[i], count_u, count_uv,
                                           0.8, 0.8);
                    adjust_params(g_trees[i], real_type(1) / g_in.total_datasets);
                }
            else if (g_in.start_trees[i] == "starz")
                {
                    vector<unsigned> parent(g_in.num_variables, 0);
                    g_trees[i].reset(parent, g_root);
                    set_params_from_counts(g_trees[i], count_u, count_uv,
                                           0.8, 0.8);
                    adjust_params(g_trees[i], real_type(1) / g_in.total_datasets);
                }
            else
                {
                    try {
                        Hot &hot = g_trees[i];
                        errno = 0;
                        FILE *infile = fopen(g_in.start_trees[i].c_str(), "r");
                        if (errno)
                            throw Exception("Could not open start tree file");

                        NHNode *root = NH_read_tree(infile, stderr, 0);
                        if (root == 0)
                            throw Exception("Could not read start tree file");

                        create_hot_from_NHNode(root, hot);
                        NH_destroy_tree(root);
                        fclose(infile);

                        adjust_params(hot, real_type(1) / g_in.total_datasets);
                    }
                    catch (Exception &e) {
                        cerr << PROG_NAME << ": Bad option given for --start-tree\n"
                             << PROG_NAME << ": " << e.what() << "\n"
                             << "Try '" << PROG_NAME << " --help' for more information.\n";
                        exit(EXIT_FAILURE);
                    }

                }
        }

    g_gamma.resize(g_in.num_trees);
    g_lambda.resize(g_in.num_trees);
    if (g_in.num_trees == 1)
        {
            g_gamma[0].resize(g_in.num_datasets, real_type(1.0));
            g_lambda[0] = real_type(1.0);
            Hot &hot = g_trees[0];

            find_tree<real_type>(0);
//             find_tree<real_type>(g_in.column_masks[0],
//                                  g_in.num_datasets,
//                                  g_in.num_variables,
//                                  0,
//                                  g_gamma[0],
//                                  g_in.max_tree_iterations,
//                                  g_in.tree_cutoff,
//                                  g_in.p_z1_zp0,
//                                  g_in.p_x1_z0,
//                                  hot);

            cout.precision(4);
            print_tree(cout, hot, 0);
            cout << ";\n";
            return EXIT_SUCCESS;
        }
    else if (g_options["hard-em"].as<bool>())
        {
            unsigned iterations = 1;
            real_type L = 0;
            real_type old_L = 0;
            /* p_X_T[x][i] == Pr[ data[x] | T_i ]*/
            multi_array<real_type, 2> p_X_T(boost::extents[g_in.num_datasets][g_in.num_trees]);
            for (unsigned i = 0; i < g_gamma.size(); ++i)
                {
                    g_gamma[i].resize(g_in.num_datasets);
                }

            /* Compute p_X_T for the all trees. */
            for (unsigned x = 0; x < g_in.num_datasets; ++x)
                {
                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            p_X_T[x][i] = compute_p_X(g_trees[i], g_in.data, x);
                        }
                }

            for (;;)
                {
                    /* Assign each datapoint to exactly one tree. This is
                       implicitly done using the gamma values sent to the
                       find_tree functions. */
                    vector<real_type> probs;

                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            unsigned max_idx = 0;
                            for (unsigned i = 1; i < g_in.num_trees; ++i)
                                {
                                    if (p_X_T[x][i] > p_X_T[x][max_idx])
                                        max_idx = i;
                                }
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                g_gamma[i][x] = real_type(0);
                            //g_gamma[max_idx][x] = p_X_T[x][max_idx];
                            if (p_X_T[x][max_idx] > real_type(1e-4))
                                g_gamma[max_idx][x] = real_type(1);

                            probs.push_back(p_X_T[x][max_idx]);
                        }
                    sort(probs.begin(), probs.end());
                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            cerr << "---\t" << probs[x] << "\n";
                        }

                    /* Lambda now holds the percentage of the assigned
                       datapoints to each tree */
                    cerr << "---\t";
                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            int num_assigned = 0;
                            for (unsigned x = 0; x < g_in.num_datasets; ++x)
                                {
                                    if (g_gamma[i][x] > real_type(0))
                                        num_assigned += g_in.input_freq[x];
                                }
                            g_lambda[i] = real_type(num_assigned) / g_in.total_datasets;
                            cerr << g_lambda[i] << "\t";
                        }


                    cerr << "\n--- " ;
                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            int sum = 0;
                            for (unsigned x = 0; x < g_in.num_datasets; ++x)
                                {
                                    if (g_gamma[i][x] > real_type(0))
                                        sum++;
                                }
                            cerr << sum << "\t";
                        }
                    cerr << "\n";

                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            find_tree<real_type>(i);
//                             find_tree<real_type>(g_in.column_masks[i],
//                                                  g_in.num_datasets,
//                                                  g_in.num_variables,
//                                                  0,
//                                                  g_gamma[i],
//                                                  g_in.max_tree_iterations,
//                                                  g_in.tree_cutoff,
//                                                  g_in.p_z1_zp0,
//                                                  g_in.p_x1_z0,
//                                                  g_trees[i]);
                        }

                    /*
                     * Compute p_X_T for the new trees.
                     */
                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    p_X_T[x][i] = compute_p_X(g_trees[i], g_in.data, x);
                                }
                        }
                    /*
                     * Compute the total likelihood.
                     */
                    L = 1;
                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            real_type p = 0;
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    p += g_gamma[i][x] * p_X_T[x][i];
                                }
                            if (p == real_type(0))
                                continue;
                            for (unsigned i = 0; i < g_in.input_freq[x]; ++i)
                                {
                                    L *= p;
                                }
                        }

                    clog << "--- old mixture likelihood: " << old_L << "\t"
                         << "new mixture likelihood: " << L << " ---\n\n";

                    /* Break the loop if cutoff reached. */
                    if (iterations >= g_in.max_mixture_iterations)
                        break;

                    real_type diff;
                    if (L < old_L)
                        {
                            clog << "--- warning: mixture likelihood decreased! ---\n";
                            diff = old_L - L;
                        }
                    else
                        {
                            diff = L - old_L;
                        }

                    if (diff/old_L < real_type(g_in.mixture_cutoff))
                        break;

                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            Hot &hot = g_trees[i];
                            clog << "--- " << g_lambda[i] << "\n";
                            clog << "--- ";
                            for (unsigned k = 0; k < hot.size(); ++k)
                                if (k != g_root)
                                    std::clog << setw(3) << k;
                            std::clog << "\n--- ";
                            for (unsigned k = 0; k < hot.size(); ++k)
                                if (k != g_root)
                                    std::clog << setw(3) << hot.parent(k);
                            std::clog << "\n";
                        }


                    /* Save old likelihood. */
                    old_L = L;
                    iterations++;
                }

            for (unsigned i = 0; i < g_in.num_trees; ++i)
                {
                    cout.precision(4);
                    cout << g_lambda[i] << "\n";
                    print_tree(cout, g_trees[i], 0);
                    cout << ";\n";
                }


        }
    else
        {
            /* Initialize the prior probabilities of the trees. */
            for (unsigned i = 0; i < g_lambda.size(); ++i)
                {
                    g_lambda[i] = real_type(g_in.start_lambdas[i]);
                }

            for (unsigned i = 0; i < g_gamma.size(); ++i)
                {
                    g_gamma[i].resize(g_in.num_datasets);
                }

            /* p_X_T[x][i] == Pr[ data[x] | T_i ]*/
            multi_array<real_type, 2> p_X_T(boost::extents[g_in.num_datasets][g_in.num_trees]);

            unsigned iterations = 1;
            real_type L = 0;
            real_type old_L = 0;
            for (;;)
                {
                    /*
                     * Compute gamma
                     */
                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    p_X_T[x][i] = compute_p_X(g_trees[i], g_in.data, x);
                                }

                            real_type sum = 0;
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    sum += g_lambda[i] * p_X_T[x][i];
                                }

                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    g_gamma[i][x] = g_lambda[i] * p_X_T[x][i] / sum;
                                }
                        }

                    /*
                     * Get the new lambdas
                     */
                    vector<real_type> gamma_sum(g_in.num_trees);
                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            gamma_sum[i] = 0;
                            for (unsigned x = 0; x < g_in.num_datasets; ++x)
                                {
                                    gamma_sum[i] += g_gamma[i][x] * g_in.input_freq[x];
                                }
                        }
                    real_type sum = 0;
                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            sum += gamma_sum[i];
                        }

                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            g_lambda[i] = gamma_sum[i] / sum;
                        }

                    /*
                     * Get the new trees.
                     */
                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            find_tree<real_type>(i);
//                             find_tree<real_type>(g_in.column_masks[i],
//                                                  g_in.num_datasets,
//                                                  g_in.num_variables,
//                                                  0,
//                                                  g_gamma[i],
//                                                  g_in.max_tree_iterations,
//                                                  g_in.tree_cutoff,
//                                                  g_in.p_z1_zp0,
//                                                  g_in.p_x1_z0,
//                                                  g_trees[i]);
                        }

                    /*
                     * Compute p_X_T for the new trees.
                     */
                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    p_X_T[x][i] = compute_p_X(g_trees[i], g_in.data, x);
                                }
                        }
                    /*
                     * Compute the total likelihood.
                     */
                    L = 1;
                    for (unsigned x = 0; x < g_in.num_datasets; ++x)
                        {
                            real_type p = 0;
                            for (unsigned i = 0; i < g_in.num_trees; ++i)
                                {
                                    p += g_lambda[i] * p_X_T[x][i];
                                }
                            for (unsigned i = 0; i < g_in.input_freq[x]; ++i)
                                {
                                    L *= p;
                                }
                        }

                    clog << "--- old mixture likelihood: " << old_L << "\t"
                         << "new mixture likelihood: " << L << " ---\n\n";

                    /* Break the loop if cutoff reached. */
                    if (iterations >= g_in.max_mixture_iterations)
                        break;

                    real_type diff;
                    if (L < old_L)
                        {
                            clog << "--- warning: mixture likelihood decreased! ---\n";
                            diff = old_L - L;
                        }
                    else
                        {
                            diff = L - old_L;
                        }

                    if (diff/old_L < real_type(g_in.mixture_cutoff))
                        break;

                    for (unsigned i = 0; i < g_in.num_trees; ++i)
                        {
                            Hot &hot = g_trees[i];
                            clog << "--- " << g_lambda[i] << "\n";
                            clog << "--- ";
                            for (unsigned k = 0; k < hot.size(); ++k)
                                if (k != g_root)
                                    std::clog << setw(3) << k;
                            std::clog << "\n--- ";
                            for (unsigned k = 0; k < hot.size(); ++k)
                                if (k != g_root)
                                    std::clog << setw(3) << hot.parent(k);
                            std::clog << "\n";
                        }


                    /* Save old likelihood. */
                    old_L = L;
                    iterations++;
                }

            for (unsigned i = 0; i < g_in.num_trees; ++i)
                {
                    cout.precision(4);
                    cout << g_lambda[i] << "\n";
                    print_tree(cout, g_trees[i], 0);
                    cout << ";\n";
                }
        }
    return EXIT_SUCCESS;
}

//*****************************************************************************
//                     Helper function declarations
//*****************************************************************************



//*****************************************************************************
//                   Function and member definitions.
//*****************************************************************************

void
set_params_from_counts(Hot &hot,
                       const multi_array<unsigned, 2> &count_u,
                       const multi_array<unsigned, 4> &count_uv,
                       real_type p_x0_z0,
                       real_type p_x1_z1)
{
    /* First set the root's parameters. */

    Param<real_type> &root_param = hot.data(hot.root());
    root_param.p_z_zp[0][0] = real_type(0);
    root_param.p_z_zp[1][0] = real_type(0);
    root_param.p_z_zp[0][1] = real_type(0);
    root_param.p_z_zp[1][1] = real_type(1);
    root_param.p_x_z[0][0] = real_type(1);
    root_param.p_x_z[0][1] = real_type(0);
    root_param.p_x_z[1][0] = real_type(0);
    root_param.p_x_z[1][1] = real_type(1);

    /* Set the z-parameters of the root's children. */
    const vector<Hot::vid_t> &root_children = hot.children(hot.root());
    for (unsigned i = 0; i < root_children.size(); ++i)
        {
            unsigned u = root_children[i];
            Param<real_type> &param = hot.data(u);
            param.p_z_zp[0][0] = real_type(0);
            param.p_z_zp[1][0] = real_type(0);
            param.p_z_zp[1][1] =
                real_type(count_u[u][1]) / real_type(g_in.total_datasets);
            param.p_z_zp[0][1] =
                real_type(count_u[u][0]) / real_type(g_in.total_datasets);
        }

    /* Set the z-parameters of the rest of the vertices. */
    for (unsigned u = 0; u < hot.size(); ++u)
        {
            if (u == hot.root() || hot.parent(u) == hot.root())
                continue;

            Param<real_type> &param = hot.data(u);
            unsigned parent = hot.parent(u);

            for (unsigned i = 0; i < 2; ++i)
                for (unsigned j = 0; j < 2; ++j)
                    {
                        if (count_u[parent][j] == 0)
                            {
                                param.p_z_zp[i][j] = real_type(0.5);
                            }
                        else
                            {
                                param.p_z_zp[i][j] =
                                    real_type(count_uv[parent][u][j][i]) /
                                    real_type(count_u[parent][j]);
                            }
                    }
        }

    /* Set the x-parameters. */
    for (unsigned u = 0; u < hot.size(); ++u)
        {
            if (u == hot.root())
                continue;

            Param<real_type> &param = hot.data(u);

            param.p_x_z[0][0] = p_x0_z0;
            param.p_x_z[1][1] = p_x1_z1;
            param.p_x_z[0][1] = real_type(1) - param.p_x_z[1][1];
            param.p_x_z[1][0] = real_type(1) - param.p_x_z[0][0];
        }
}


void
adjust_params(Hot &hot, real_type zero_replacement)
{
    /* Ensure that the root has the special parameter values. */
    Param<real_type> &root_param = hot.data(hot.root());
    root_param.p_z_zp[0][0] = real_type(0);
    root_param.p_z_zp[1][0] = real_type(0);
    root_param.p_z_zp[0][1] = real_type(0);
    root_param.p_z_zp[1][1] = real_type(1);
    root_param.p_x_z[0][0] = real_type(1);
    root_param.p_x_z[0][1] = real_type(0);
    root_param.p_x_z[1][0] = real_type(0);
    root_param.p_x_z[1][1] = real_type(1);

    /* ensure that no probability is zero. */
    for (unsigned u = 0; u < hot.size(); ++u)
        {
            if (u == hot.root())
                continue;

            Param<real_type> &param = hot.data(u);

            for (unsigned i = 0; i < 2; ++i)
                for (unsigned j = 0; j < 2; ++j)
                    {
                        if (param.p_z_zp[i][j] < zero_replacement)
                            {
                                param.p_z_zp[i][j] = zero_replacement;
                                param.p_z_zp[1-i][j] =
                                    real_type(1) - zero_replacement;
                            }
                        if (param.p_x_z[i][j] < zero_replacement)
                            {
                                param.p_x_z[i][j] = zero_replacement;
                                param.p_x_z[1-i][j] =
                                    real_type(1) - zero_replacement;
                            }
                    }
        }

    /* ensure that the root's children have the required special values. */
    const vector<Hot::vid_t> &root_children = hot.children(hot.root());
    for (unsigned i = 0; i < root_children.size(); ++i)
        {
            unsigned u = root_children[i];
            Param<real_type> &param = hot.data(u);
            param.p_z_zp[0][0] = real_type(0);
            param.p_z_zp[1][0] = real_type(0);
        }
}


void
extend_parent(vector<unsigned> &parent,
              dynamic_bitset<> mask,
              unsigned root)
{
    unsigned small_size = parent.size();

    vector<unsigned> mapping(small_size);
    unsigned nr = 0;
    for (dynamic_bitset<>::size_type i = mask.find_first(); i != dynamic_bitset<>::npos;
         i = mask.find_next(i))
    {
        mapping[nr++] = i;
    }

    vector<unsigned> new_parent(mask.size(), root);
    for (unsigned i = 0; i < small_size; ++i)
        {
            new_parent[mapping[i]] = mapping[parent[i]];
        }

    parent.swap(new_parent);
}
