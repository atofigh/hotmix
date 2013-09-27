#include <vector>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <algorithm>
#include <stack>

template<class E, class A> inline
void
Assert(A assertion)
{
    if (!assertion) throw E();
}

/*
 * Static variable definitions.
 */

template<typename D>
const typename Rooted_tree<D>::vid_t Rooted_tree<D>::NONE;

/*
 * Definitions of exception classes used in Rooted_tree
 */

template<typename D>
const char *
Rooted_tree<D>::Invalid_id::what() const throw()
{
    return "Invalid vertex id passed to Rooted_tree class.";
}

template<typename D>
const char *
Rooted_tree<D>::Invalid_parent::what() const throw()
{
    return "Invalid parent specification passed to Rooted_tree.";
}

/*
 * Constructor and functions definitions.
 */

template<typename D>
Rooted_tree<D>::Rooted_tree()
    : root_(NONE), postorder_begin_(NONE)
{    
}

template<typename D>
Rooted_tree<D>::Rooted_tree(const std::vector<vid_t> &parent, vid_t root)
{
    reset(parent, root);
}

template<typename D>
void
Rooted_tree<D>::reset(const std::vector<vid_t> &parent, vid_t root)
{
    Assert<Invalid_parent>(!parent.empty());

    root_ = root;
    vertices_.clear();
    vertices_.resize(parent.size());
    preorder_next_.clear();
    postorder_next_.clear();
    E_.clear();
    L_.clear();
    Ref_.clear();
    
    /* Fill in the parent and children information. */
    vertices_[root].parent = NONE;
    for (unsigned i = 0; i < parent.size(); ++i)
        {
            if (i == root_)
                continue;

            Assert<Invalid_parent>(parent[i] < parent.size());
            vertices_[i].parent = parent[i];
            vertices_[ parent[i] ].children.push_back(i);
        }

    /* Build lca datastructures. */
    build_lca_();
    
    /* create the preorder and postorder traversal data structures. */
    std::vector<vid_t> preorder;
    std::stack<vid_t> s;
    s.push(root_);
    while (!s.empty())
        {
            vid_t cur = s.top(); s.pop();
            preorder.push_back(cur);
            BOOST_FOREACH (vid_t child, children(cur))
                {
                    s.push(child);
                }
        }
    preorder_next_.resize(size());
    postorder_next_.resize(size());

    postorder_begin_ = preorder.back();
    preorder_next_[preorder[0]] = preorder[1];
    preorder_next_[preorder.back()] = NONE;
    postorder_next_[preorder[0]] = NONE;
    postorder_next_[preorder.back()] = preorder[preorder.size()-2];
    for (unsigned i = 1; i < preorder.size() - 1; ++i)
        {
            preorder_next_[preorder[i]] = preorder[i+1];
            postorder_next_[preorder[i]] = preorder[i-1];
        }
}

     


template<typename D>
void
Rooted_tree<D>::create_EL_(vid_t cur_vertex, unsigned cur_level)
{
    E_.push_back(cur_vertex);
    L_.push_back(cur_level);
            
    if (this->is_leaf(cur_vertex))
        return;
    
    BOOST_FOREACH (unsigned child, vertices_[cur_vertex].children)
        {
            create_EL_(child, cur_level + 1);
            E_.push_back(cur_vertex);
            L_.push_back(cur_level);
        }
}


template<typename D>
template<typename T>
unsigned
Rooted_tree<D>::msb_(T v)
{
    static const char LogTable256[] = 
        {
            0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
    
    unsigned bits = std::numeric_limits<T>::digits;
    unsigned c = 0; // c will be lg(v) truncated
    while (bits > 8)
        {
            bits /= 2;
            if (v >> bits)
                {
                    c += bits;
                    v >>= bits;
                }
        }
    c += LogTable256[v];
    return c;
}


template<typename D>
void
Rooted_tree<D>::build_lca_()
{
    const unsigned N = this->size();

    create_EL_(this->root(), 0);
    Assert<Invalid_parent>(L_.size() == 2*N - 1);

    M_.resize(boost::extents[L_.size()][msb_(L_.size()) + 1]);
    Ref_.resize(this->size(), NONE);
    for (unsigned i = 0; i < E_.size(); ++i)
        {
            if (Ref_[E_[i]] > N)
                Ref_[E_[i]] = i;
        }

    unsigned rows = L_.size();
    unsigned cols = msb_(L_.size()) + 1;

    for (unsigned i = 0; i < rows; ++i)
        M_[i][0] = i;
    
    for (unsigned j = 1; j < cols; ++j)
        {
            unsigned stride = 1u << j;
            for (unsigned i = 0; i < rows; ++i)
                {
                    if (i + stride/2 >= rows)
                        {
                            M_[i][j] = M_[i][j-1];
                            continue;
                        }
                    unsigned idx1 = M_[i][j-1];
                    unsigned idx2 = M_[i + stride/2][j-1];
                    
                    M_[i][j] = L_[idx2] < L_[idx1] ? idx2 : idx1;
                }
        }
}


template<typename D>
typename Rooted_tree<D>::size_type
Rooted_tree<D>::size() const
{
    return vertices_.size();
}


template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::root() const
{
    return root_;
}


template<typename D>
const std::vector< typename Rooted_tree<D>::vid_t > &
Rooted_tree<D>::children(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices_[v].children;
}

template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::parent(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices_[v].parent;
}

template<typename D>
const typename Rooted_tree<D>::data_type &
Rooted_tree<D>::data(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices_[v].data;
}


template<typename D>
typename Rooted_tree<D>::data_type &
Rooted_tree<D>::data(vid_t v)
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices_[v].data;
}


template<typename D>
bool
Rooted_tree<D>::is_leaf(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return vertices_[v].children.empty();
}

template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::preorder_begin() const
{
    return root();
}

template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::preorder_next(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return preorder_next_[v];
}
 
template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::postorder_begin() const
{
    return postorder_begin_;
}

template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::postorder_next(vid_t v) const
{
    Assert<Invalid_id>(v >= 0 && v < size());

    return postorder_next_[v];
}

/*
 * The algorithm is from the article 'Lowest common ancestors in trees
 * and directed acyclic graphs' by M Bender et al, 2005. The tables
 * and arrays have the same names as in the article, except for the
 * R-array which here is renamed to Ref.
 */
template<typename D>
typename Rooted_tree<D>::vid_t
Rooted_tree<D>::lca(vid_t v1, vid_t v2) const
{
    Assert<Invalid_id>(v1 >= 0 && v1 < size());
    Assert<Invalid_id>(v2 >= 0 && v2 < size());

    /* Get the representatives (i.e, indexes into L) of v1 and v2. */
    unsigned r1 = Ref_[v1];
    unsigned r2 = Ref_[v2];

    /* Make sure that r2 is the bigger one, so that the range
       of indices is [r1, r2]. */
    if (r1 > r2)
        std::swap(r1, r2);

    unsigned k = msb_(r2 - r1 + 1);
    unsigned idx1 = M_[r1][k];
    unsigned idx2 = M_[r2-(1u<<k)+1][k]; /* 1u<<k == 2^k */
    if (L_[idx2] < L_[idx1])
        idx1 = idx2;

    return E_[idx1];
}
    

template<typename D>
bool
Rooted_tree<D>::descendant(vid_t v1, vid_t v2) const
{
    /* Is v1 descendant of v2? */

    return lca(v1, v2) == v2;
}
