#ifndef COMMON_HPP
#define COMMON_HPP

extern "C" {
#include <NHparser/NHparser.h>
}
#include "rooted-tree.hpp"
#include <mpfrcpp/mpfrcpp.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <exception>
#include <vector>
#include <string>

typedef mpfr::Real MP_double;

/*
 * Class Param<Real>
 *
 * An object of type Param keeps track of the parameters of an edge in
 * a HOT.
 */
template<typename Real>
struct Param {
    Real p_x_z[2][2];   /* p_x_z[i][j]  = Pr[X(u)=i | Z(u)=j, HOT] */
    Real p_z_zp[2][2];  /* p_z_zp[i][j] = Pr[Z(u)=i | Z(p(u))=j, HOT] */
};

/*
 * class Exception
 *
 * Exception class from which to derive all the exceptions.
 */

class Exception : std::exception {
public:
    Exception(std::string message) : message_(message) {}
    ~Exception() throw() {}
    const char *what() const throw() {return message_.c_str();}
private:
    std::string message_;
};

/*
 * Class Bad_hot
 *
 * An exception class that might be thrown by create_hot_from_NHNode().
 */
struct Bad_hot : public Exception {
    Bad_hot(std::string message) : Exception(message) {}
};

/*
 * Class Bad_mixture_spec
 *
 * An exception class that might be thrown by read_mixture_spec().
 */
struct Bad_mixture : public Exception {
    Bad_mixture(std::string message) : Exception(message) {}
};

/*
 * init_rand<>()
 *
 * Seeds the random number generator 'engine'. Assumes engine has a
 * function named 'seed', taking an unsigned.
 */
template <class RandomNumberGenerator>
void init_rand(RandomNumberGenerator &engine);

template <class RandomNumberGenerator>
void init_rand(RandomNumberGenerator &engine, unsigned seed);

/*
 * log(Real)
 *
 * We need to define a log operation for all the Real-types we use.
 */
inline MP_double log(MP_double d);

/*
 * gen_random_tree<>()
 *
 * Creates a random tree without setting the parameters. The result is
 * a parent vector that may be passed to a Rooted_tree.
 */
template<typename Real>
void
gen_random_tree(unsigned N,
                unsigned root,
                std::vector<unsigned> &parent);

/*
 * print_tree<>()
 *
 * prints a HOT to the output stream 'out'. 'u' sepcifies the root of
 * the (sub)tree that should be printed which should be set to the
 * root vertex of the tree to print the entire tree. Note that no
 * semicolon or newline is printed to 'out'. Every vertex has its
 * number as label followed by a comment containing all eight
 * probabilities. The order of the probabilites is:
 *
 * Pr[Z(u)=0|Z(p(u))=0] Pr[Z(u)=0|Z(p(u))=1] Pr[Z(u)=1|Z(p(u))=0] Pr[Z(u)=1|Z(p(u))=1]
 * Pr[X(u)=0|Z(u)=0] Pr[X(u)=0|Z(u)=1] Pr[X(u)=1|Z(u)=0] Pr[X(u)=1|Z(u)=1]
 *
 */
template<typename Real>
void
print_tree(std::ostream &out, Rooted_tree< Param<Real> > &tree, unsigned u);

/*
 * create_hot_from_NHNode<>()
 *
 * Given a tree with HOT-annotations, returns the corresponding HOT in
 * 'hot'. Might throw an instance of Bad_hot.
 */
template<typename Real>
void
create_hot_from_NHNode(NHNode *root,
                       Rooted_tree< Param<Real> > &hot);

/*
 * read_mixture()
 *
 * Reads a mixture and retuns the result in 'hots' and 'probs'. Throws
 * instances of Bad_mixture if anything goes wrong.
 */
template<typename Real>
void
read_mixture(std::string mixture_filename,
             std::vector< Rooted_tree< Param<Real> > > &hots,
             std::vector<Real> &probs);

/*
 * read_input_lines()
 *
 * Reads the input stream 'in' line by line. Each line is stripped of
 * any and all whitespace. Those (stripped) lines that contain only
 * the characters '1' and '0' are inserted into the vector
 * 'input_lines'.
 */
void
read_input_lines(std::istream &in,
                 std::vector<std::string> &input_lines);

/*
 * compute_p_X<>()
 *
 * computes the probability of data[dataset] given a hot.
 */
template<typename Real, typename Array_2d>
Real
compute_p_X(const Rooted_tree< Param<Real> > &hot,
            const Array_2d &data,
            unsigned dataset);

/*
 * operator<<()
 *
 * Used to print Log_doubles.
 */
// inline
// std::ostream &
// operator<<(std::ostream &out, const Log_double &ld);

/*
 * Global random number generators:
 *
 * 'g_rng_ui': used to get a uniformly distributed random
 * unsigned. g_rng_ui(n) gives a random number in the range [0,n).
 *
 * 'g_rng_d': used to get a uniformly distributed random double in the
 * range [0,1).
 *
 * The random number engine 'g_generator' is seeded in main.
 */
extern boost::mt19937    g_generator;
extern boost::random_number_generator<boost::mt19937, unsigned> g_rng_ui;
extern boost::variate_generator<boost::mt19937&, boost::uniform_real<> > g_rng_d;

/*
 * Default values for options passed to the different programs in the
 * package.
 */
const double g_p_z1_zp0 = 1e-3;
const double g_p_x1_z0 = 1e-3;
const double g_min_z_prob = 1e-3;
const double g_min_x_prob = 1e-3;
const double g_min_mix_prob = 1e-3;
const double g_sum_tolerance = 1e-3;
const double g_min_p_z1_zp1 = 0.5;
const double g_min_p_x1_z1 = 0.5;

#include "common-impl.hpp"

#endif /* COMMON_HPP */
