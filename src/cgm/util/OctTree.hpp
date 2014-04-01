//-------------------------------------------------------------------------
// Filename      : NodeSearch.hpp
//
// Purpose       : Class for O(ln n) search for nodes within a specified
//                 tolerance of a specified position.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------

#ifndef OCT_TREE_HPP
#define OCT_TREE_HPP

#include "CubitVector.hpp"
#include "DLIList.hpp"

template <class X, class E> class OctTreeCell;
template <class X, class E> struct OctTreeEntry;

template <class X> class DefaultPointQuery{
  public: static inline CubitVector coordinates(const X* x)
    { return x->coordinates(); }
};

template < class X, class E = DefaultPointQuery<X> >
class OctTree
{

  public:
  
    OctTree( DLIList<X*>& nodes, double tolerance,
                int min_nodes_per_box = -1,
                double min_box_dimension = -1.0 );
      //- Constructor
      //- nodes     - list of nodes to build oct-tree for
      //- tolerance - tolerance for comparison of node positions
      //- min_nodes_per_box - don't subdivide boxes in tree with 
      //-                     this many or fewer nodes.
      //- min_box_dimension - don't subdivide boxes in tree if
      //-                     all three dimensions are less than this
    
    ~OctTree();
      // - Destructor
    
    void nodes_near( const CubitVector& position,
                     DLIList<X*>& result_list );
      //- Get the list of all nodes in the tree for which 
      //- the distance between the passed position and the 
      //- node is less than or equal to the tolerance 
      //- value passed to the constructor.
    
    void tree_size( DLIList<int>& count_at_each_level,
                    DLIList<int>& leaves_at_each_level );
      //- Get the size of the tree at each level.
    
  private:
  
    double tolerance_;
    double min_box_size_;
    int min_nodes_;
      // values passed to constructor
  
    OctTreeCell<X,E>* root_;
      // root of oct-tree
    
    CubitVector min_;
    CubitVector max_;
      // bounding box of entire oct-tree (and size of root box)
    
    OctTreeEntry<X,E>* node_memory_;
      // memory allocated for linked-list entries
      // this is allocated by the constructor, passed to
      // the root node of the oct-tree for its internal use
      // and released by the destructor.  All use of this
      // memory is internal to the NodeSearchBox class.

    void split_node( OctTreeCell<X,E>* node,
                     const CubitVector& min,
                     const CubitVector& max );
      //- Recursive function to subdivide octree nodes
  
    OctTreeCell<X,E>* allocate_8();
      //- Allocate a block of 8 oct-tree nodes for use as
      //- children of some current node.
    
    DLIList<OctTreeCell<X,E>*> mem_pages_;
    OctTreeCell<X,E>* curr_page_end_;
      //- This data is the memory pool from which allocate_8()
      //- allocates oct-tree nodes.
  
    void push_search_box( int quadrant_flags );
      //- Push the child box of 'current_box' indicated
      //- by quadrant_flags onto box_list stack or onto
      //- leaf_list.
    
    void pop_to_current();
      //- Pop one box from box_list and store it in 
      //- current_* variables.  Also updates 'center'
    
    OctTreeCell<X,E>* box_list[64];  // stack of none-leaf boxes to be searched
    CubitVector box_min[64]; // minimum of boxes in box_list
    CubitVector box_max[64]; // maximum of boxes in box_list
    int box_count;          // count of boxes in box_list
    DLIList<X*> search_results; // results of a search
      // State information for tree search
    
    CubitVector current_min;
    CubitVector current_max;
    CubitVector center;
    OctTreeCell<X,E>* current_box;
      // More state information for tree search
      // The 'current' box during a search, and related info.
      
      
    void tree_size( OctTreeCell<X,E>* box,
                    DLIList<int>& boxes,
                    DLIList<int>& leaves );
                    
};

#ifdef TEMPLATE_DEFS_INCLUDED
#  define INCLUDED_FROM_OCT_TREE_HPP
#  include "OctTree.cpp"
#  undef INCLUDED_FROM_OCT_TREE_HPP
#endif

#endif

