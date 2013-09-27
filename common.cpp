#include "common.hpp"
#include <boost/foreach.hpp>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <string>

using namespace std;

boost::mt19937           g_generator;
boost::random_number_generator<boost::mt19937, unsigned>
                         g_rng_ui(g_generator);
boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
                         g_rng_d(g_generator, boost::uniform_real<>());


void
read_input_lines(istream &in, vector<string> &input_lines)
{
    /* Read the stream line by line. */
    string line;
    while (getline(in, line))
        {
            /* strip the line from all whitespace */
            string stripped;
            boost::tokenizer<> tok(line);
            BOOST_FOREACH (string token, tok)
                {
                    stripped += token;
                }

            if (stripped.size() == 0)
                continue;

            /* keep stripped line if it contains only zeros and ones. */
            if (stripped.find_first_not_of("01") == string::npos)
                input_lines.push_back(stripped);
        }
}
