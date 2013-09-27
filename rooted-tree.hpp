#ifndef ROOTED_TREE_HPP
#define ROOTED_TREE_HPP

#include <exception>
#include <vector>
#include <boost/multi_array.hpp>

/*
 * [ type Rooted_tree ]: The class Rooted_tree implements the rooted
 * tree concept. Associated with each vertex is a field of type
 * 'VertexDataType'.
 */

template <typename VertexDataType>
class Rooted_tree {
public:
    /*
     * [ type vid_t ]: An unsigned integral type used for vertex
     * ids. A constant value 'NONE' is reserved to signal a
     * non-existant vertex id and is guaranteed not to be allocated
     * for a vertex.
     */
    typedef unsigned                            vid_t;

    typedef std::vector<vid_t>::size_type       size_type;
    typedef VertexDataType                      data_type;

    /*
     * We cannot use numeric_limits<vid_t>::max() to initialize NONE.
     * It is illegal to call a function here, and when initializing it
     * outside, we run into a bug in gcc.
     */
    static const vid_t                          NONE = -1;

    /*
     * [ Exception classes ]: Exceptions that may be thrown by the
     * functions in this class. All derived from std::exception.
     */
    struct Invalid_id 
        : public std::exception {const char *what() const throw();};
    struct Invalid_parent : 
        public std::exception {const char *what() const throw();};
    
    /*
     * [ Constructors ]: Constructs a tree with parent.size()
     * vertices, where the parent of vertex u is given by
     * parent[u]. The root vertex is specified by 'root'. Also builds
     * datastructures used for computing lca. Therefore, the
     * construction takes time O(n log n).
     *
     * VertexDataType must be default constructible.
     *
     * Default Copying and assignment are used.
     */
    Rooted_tree();
    Rooted_tree(const std::vector<vid_t> &parent, vid_t root);
    void reset(const std::vector<vid_t> &parent, vid_t root);

    /*
     * [ size ]: Returns the number of vertices in the tree.
     */
    size_type size() const;

    /*
     * [ root ]: Returns the root vertex id.
     */
    vid_t root() const;

    /*
     * [ children ]: returns a vector containing the children of
     * 'vertex'.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    const std::vector<vid_t> &children(vid_t vertex) const;

    /*
     * [ parent ]: returns the parent of 'vertex', or NONE if 'vertex'
     * is the root.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    vid_t parent(vid_t vertex) const;

    /*
     * [ data ]: Returns the data associated with 'vertex'.
     * 
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    const data_type &data(vid_t vertex) const;
    data_type &data(vid_t vertex);

    /*
     * [ is_leaf ]: Returns true if vertex has no children.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    bool is_leaf(vid_t vertex) const;

    /*
     * [ preorder_begin/preorder_next]: preorder_begin returns the
     * first vertex of a preorder traversal of the tree (this is
     * always the root, but was included for symmetri with the
     * postorder functions). preorder_next returns the next vertex in
     * a preorder traversal or NONE if there are no more vertices to visit.
     *
     * A preorder traversal of the tree
     *        ---4
     *       |
     *   ----2
     *   |   |
     *---0    ---3
     *   |
     *   --------1
     *
     * could be 0, 1, 2, 3, 4.
     *
     * The time it takes for these functions to return is not
     * constant, but the total contribution from these two functions
     * in a loop that visits all vertices is O(n), where n is the
     * number of vertices in the tree.
     *
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    vid_t preorder_begin() const;
    vid_t preorder_next(vid_t vertex) const;

    /*
     * [ postorder_begin/postorder_next]: postorder_begin returns the
     * first vertex of a postorder traversal of the
     * tree. postorder_next returns the next vertex in a postorder
     * traversal or NONE if there are no more vertices to visit.
     *
     * A postorder traversal of the tree
     *        ---4
     *       |
     *   ----2
     *   |   |
     *---0    ---3
     *   |
     *   --------1
     *
     * could be 4, 3, 2, 1, 0
     *
     * The comments on time complexity for the preorder-functions
     * above apply for the postorder-functions as well.
     * 
     * Precondition:
     * 1. 'vertex' is a valid vertex id.
     */
    vid_t postorder_begin() const;
    vid_t postorder_next(vid_t vertex) const;

    /*
     * [ lca ]: Returns the the least common ancestor of vertices v1
     * and v2. The first time this function is called, some
     * precomputation is done that takes time O(n log n). Thereafter
     * the lca is computed in constant time. Note, however, that after
     * each change to the tree (e.g. after a call to add_children),
     * the precomputation has to be redone!
     *
     * The algorithm is from the article 'Lowest common ancestors in
     * trees and directed acyclic graphs' by M Bender et al, 2005.
     *
     * Precondition:
     * 1. 'v1' and 'v2' are valid vertex ids.
     */
    vid_t lca(vid_t v1, vid_t v2) const;

    /*
     * [ descendant ]: Returns true if 'v1' is a descendant of
     * 'v2'. This function is simply a wrapper around lca, since 'v1'
     * is a descendant of 'v2' if and only if lca(v1,v2) = v2.
     */
    bool descendant(vid_t v1, vid_t v2) const;


private:
    struct vertex_t {
        vid_t parent;
        std::vector<vid_t> children;
        data_type data;
    };

    /* Functions for building lca. */
    void build_lca_();
    void create_EL_(vid_t cur_vertex, unsigned level);
    template<typename T> static unsigned msb_(T v);
    
    vid_t root_;
    std::vector<vertex_t> vertices_;
    std::vector<vid_t> preorder_next_;
    std::vector<vid_t> postorder_next_;
    vid_t postorder_begin_;
    std::vector<vid_t> E_; /* Euler-path for lca-compuatation. */
    std::vector<unsigned> L_; /* Level array corresponding to E. */
    std::vector<std::vector<vid_t>::size_type> Ref_;
                             /* Representative array for lca-computation */
    boost::multi_array<std::vector<vid_t>::size_type, 2> M_; 
                                     /* M-matrix for the RMQ-algorithm. */
};


#include "rooted-tree-impl.hpp"


#endif /* not ROOTED_TREE_HPP. */
