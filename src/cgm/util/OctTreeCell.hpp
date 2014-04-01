//-------------------------------------------------------------------------
// Filename      : OctTreeCell.hpp
//
// Purpose       : Oct-tree node used by OctTree
//
// Special Notes : OctTree handles all memory allocation for this
//                 class.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------

#ifndef OCT_TREE_CELL_HPP
#define OCT_TREE_CELL_HPP

template <class X, class E> class OctTreeCell;
class CubitVector;
template <class X> class DLIList;

// Linked-list node used to hold CubitNodes in each box.
template <class X, class E> struct OctTreeEntry
{
  X* node;
  OctTreeEntry<X,E>* next;
};

// A node in the oct-tree.
template <class X, class E> class OctTreeCell
{
  private:
    
    /*
    union {
      OctTreeEntry<X,E>* head_; 	// head of linked-list of CubitNodes
      OctTreeCell<X,E>* children_; // array of 8 child boxes
    };
    */
    // HP-UX aCC A.03.30 chokes on this union
    void* head_or_children;
      // union of both the head pointer for the linked-list
      // of CubitNodes in a leaf node and a pointer to the
      // array of child boxes.  
      //
      // This node is a leaf node iff node_count_ is non-zero
      // or both values in the data union are null.
    
    int node_count_;
      // Number of nodes contained in a leaf node.  If zero,
      // need to check if data.children is null to determine
      // if this is a leaf node.
    
    void add_nodes( OctTreeEntry<X,E>* node );
      // Join the passed linked-list to the data list.
      // asserts if not a leaf node.
    
  public:
  
    OctTreeCell( OctTreeEntry<X,E> array[], int size );  
      // create a root node
      // calling code provides storage for linked-list
      // nodes by providing 'array'.  Do not delete 
      // 'array' until this OctTreeCell and all its
      // children are destroyed.
    
    inline OctTreeCell( ) : head_or_children(0), node_count_(0) {}
      // Construct unused child node.  This is provieded
      // for the cexternal memory allocator.
    
    bool leaf() const { return node_count_ || !head_or_children; }
      // Is this box a leaf in the tree?
    
    bool split( const CubitVector& my_center, OctTreeCell<X,E>* memory );
      // Change from leaf box to non-leaf.  Splits nodes amongst 8 new
      // child boxes using the value of 'my_center'.  May fail (return
      // false) if two few nodes are contained for a split to be sensible.
      //
      // Calling application provides memory for use in holding child
      // nodes.  memory must be freed by the caller AFTER the destruction
      // of this object.
    
    void append_nodes( DLIList<X*>& list );
      // If this node is a leaf node, append all contained CubitNodes
    
    void append_all_nodes( DLIList<X*>& list );
      // If this node is a leaf node, append all contained CubitNodes.
      // Otherwise descend oct-tree getting CubitNodes contained in
      // all child leaf nodes.

    enum {
      X_MIN = 0,
      Y_MIN = 0,
      Z_MIN = 0,
      X_MAX = 1,
      Y_MAX = 2,
      Z_MAX = 4
    };
      // Constants for selecting one of 8 child boxes.
      // Can be either bitwise-or'ed (|) or added (+).
      // MAX values take precedence over MIN values 
      // if both are specified.  (MIN is just the
      // absense of MAX).
    
    OctTreeCell<X,E>* child( int quadrant );
      // quadrant should be a bitwise OR (or sum) of 
      // the constants above.  If both a MIN and MAX 
      // value is specified for a given principal 
      // direction, the M  : children_(0),
    
    void node_locations( CubitVector& box_min, 
                         CubitVector& box_max,
                         CubitVector& average,
                         CubitVector* std_dev = 0 );
      // Get information about CubitNodes contained in
      // this box.  Method does nothing when called on
      // non-leaf boxes.
                         
    int node_count() { return node_count_; }
      // CubitNodes contained in this box.
      // Always zero for non-leaf boxes.
    
#ifdef BOYD15
    int count_all_nodes();
      // Get nodes contained in this box and all child boxes.
#endif
};

#ifdef TEMPLATE_DEFS_INCLUDED
#  define INCLUDED_FROM_OCT_TREE_CELL_HPP
#  include "OctTreeCell.cpp"
#  undef INCLUDED_FROM_OCT_TREE_CELL_HPP
#endif


#endif

