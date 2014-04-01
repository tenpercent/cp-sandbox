//-------------------------------------------------------------------------
// Filename      : OctTree.cpp
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

#if !defined(TEMPLATE_DEFS_INCLUDED) || defined(INCLUDED_FROM_OCT_TREE_HPP)

#include "OctTree.hpp"
#include "OctTreeCell.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"

const int DEFAULT_MIN_NODES_PER_BOX = 6;
  // The default value to use for the minimum number of
  // nodes in a box, if none is specified in the constructr.

const int OCT_TREE_BOX_PAGE_SIZE = 8192;
  // The size of memory blocks to allocate for holding
  // oct-tree nodes.

const int OCT_TREE_CHUNK_SIZE = 8 * sizeof(OctTreeCell<int,int>);
  // Do not change.  Size of a child array for a oct-tree node.

const int OCT_TREE_CHUNK_PER_PAGE = OCT_TREE_BOX_PAGE_SIZE / 
                                       OCT_TREE_CHUNK_SIZE;
  // Do not change.  Number of arrays of 8 oct-tree nodes 
  // (child arrays) fit in a page.
                                       
const int OCT_TREE_ALLOC_COUNT = OCT_TREE_CHUNK_PER_PAGE * 8;
  // When allocating the page as one big array of oct-tree nodes,
  // the size of the array to request.  Will be sligthtly smaller
  // than the requested page because page size probably is not a
  // multiple of OCT_TREE_CHUNK_SIZE.  That's a good thing though,
  // because then the std-c and c++ memory allocators have a little
  // extra space for internal data.

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
OctTree<X,E>::OctTree( DLIList<X*>& nodes, 
                        double tolerance,
                        int min_nodes_per_box,
                        double min_box_dimension )
{
    // Check and store constructor arguments
  tolerance_ = tolerance < GEOMETRY_RESABS ? GEOMETRY_RESABS : tolerance;
  min_nodes_ = min_nodes_per_box < 0 ? DEFAULT_MIN_NODES_PER_BOX : min_nodes_per_box;
  double tol3 = 3 * tolerance_;
  min_box_size_ = min_box_dimension < tol3 ? tol3 : min_box_dimension;
    // Note: min_box_size must be greater than twice the tolerance
    //       or the internal stack used for searching the tree will
    //       overflow (more than 8 boxes may be within the tolerance
    //       of a passed position).
  
    // set up data for memory pool of oct-tree nodes
  mem_pages_.append( new OctTreeCell<X,E>[OCT_TREE_ALLOC_COUNT] );
  curr_page_end_ = mem_pages_.get() + OCT_TREE_ALLOC_COUNT;

    // construct root node
  node_memory_ = new OctTreeEntry<X,E>[nodes.size()];
  nodes.reset();
  OctTreeEntry<X,E>* ptr = node_memory_;
  for( int i = nodes.size(); i--; )
  {
    ptr->node = nodes.get_and_step();
    ptr->next = 0;
    ptr++;
  }
  root_ = new OctTreeCell<X,E>( node_memory_, nodes.size() );
  
    // get bounding box of all nodes
  CubitVector junk;
  root_->node_locations( min_, max_, junk );
  
    // create the oct-tree
  split_node(root_, min_, max_ );
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
OctTree<X,E>::~OctTree()
{
    // Release all memory
    // Note that the oct-tree is not deleted except for the root node.
    // Only the root node was allocated directly from the heap.  All
    // other nodes were allocated from our internal memory bool by
    // the allocate_8() method.  We just release the whole pool.
  delete root_;
  delete node_memory_;
  while( mem_pages_.size() )
    delete mem_pages_.pop();
    
    // Reinitialize to catch stale pointer to this object.
  
  node_memory_ = 0;
  root_ = 0;
}

//-------------------------------------------------------------------------
// Purpose       : Find all the nodes within tolerance of the passed position.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
void OctTree<X,E>::nodes_near( const CubitVector& position, 
                             DLIList<X*>& result_list )
{
  if( position.x() - min_.x() < -tolerance_ ||
      position.y() - min_.y() < -tolerance_ ||
      position.z() - min_.z() < -tolerance_ ||
      position.x() - max_.x() >  tolerance_ ||
      position.y() - max_.y() >  tolerance_ ||
      position.z() - max_.z() >  tolerance_  )
    return;
  
    // Initialize search stack
  if( root_->leaf() )
  {
      // Done searching already -- that was easy
    root_->append_nodes( search_results );
    box_count = 0;
  }
  else
  {
    search_results.clean_out();

      // Start with root node on stack of nodes to search
    box_count = 1;
    box_list[0] = root_;

    box_min[0] = min_;
    box_max[0] = max_;
  }
  
    // While there are still boxes on the search stack
  while( box_count )
  {
      // Pop the 'current' box from the stack
    pop_to_current();
    
      // Push appropriate child box(es) to stack
    CubitVector diff = position - center;
    if( diff.x() <= tolerance_ )
    {
      if( diff.y() <= tolerance_ )
      {
        if( diff.z() <= tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MIN |
                           OctTreeCell<X,E>::Y_MIN |
                           OctTreeCell<X,E>::Z_MIN );
        }
        if( diff.z() >= -tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MIN |
                           OctTreeCell<X,E>::Y_MIN |
                           OctTreeCell<X,E>::Z_MAX );
        }
      }
      if( diff.y() >= -tolerance_ )
      {
        if( diff.z() <= tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MIN |
                           OctTreeCell<X,E>::Y_MAX |
                           OctTreeCell<X,E>::Z_MIN );
        }
        if( diff.z() >= -tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MIN |
                           OctTreeCell<X,E>::Y_MAX |
                           OctTreeCell<X,E>::Z_MAX );
        }
      }
    }      
    if( diff.x() >= -tolerance_ )
    {
      if( diff.y() <= tolerance_ )
      {
        if( diff.z() <= tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MAX |
                           OctTreeCell<X,E>::Y_MIN |
                           OctTreeCell<X,E>::Z_MIN );
        }
        if( diff.z() >= -tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MAX |
                           OctTreeCell<X,E>::Y_MIN |
                           OctTreeCell<X,E>::Z_MAX );
        }
      }
      if( diff.y() >= -tolerance_ )
      {
        if( diff.z() <= tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MAX |
                           OctTreeCell<X,E>::Y_MAX |
                           OctTreeCell<X,E>::Z_MIN );
        }
        if( diff.z() >= -tolerance_ )
        {
          push_search_box( OctTreeCell<X,E>::X_MAX |
                           OctTreeCell<X,E>::Y_MAX |
                           OctTreeCell<X,E>::Z_MAX );
        }
      }
    }      
  }
  
    // append to result list all nodes within tolerance
    // of passed position
  search_results.reset();
  const double tol_sqr = tolerance_ * tolerance_;
  for( int i = search_results.size(); i--; )
  {
    X* node = search_results.get_and_step();
    CubitVector coords = E::coordinates(node);
    double x = position.x() - coords.x();
    double y = position.y() - coords.y();
    double z = position.z() - coords.z();
    
    double dist_sqr = x * x + y * y + z * z;
    if( dist_sqr <= tol_sqr )
      result_list.append( node );
  }
}


//-------------------------------------------------------------------------
// Purpose       : Push child boxes onto search stack (or leaf node list)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
void OctTree<X,E>::push_search_box( int quadrant_flags )
{
  OctTreeCell<X,E>* box = current_box->child( quadrant_flags );
  if( box->leaf() )
  {
      // If a leaf node, get its nodes
    box->append_nodes( search_results );
  }
  else
  {
    assert( box_count < 64 );

      // Get appropriate min and max corners of child box
      
    CubitVector& oldmin = current_min;
    CubitVector& oldmax = current_max;
    CubitVector& newmin = box_min[box_count];
    CubitVector& newmax = box_max[box_count];
    
    if( quadrant_flags & OctTreeCell<X,E>::X_MAX )
    {
      newmin.x( center.x() );
      newmax.x( oldmax.x() );
    }
    else
    {
      newmin.x( oldmin.x() );
      newmax.x( center.x() );
    }
    
    if( quadrant_flags & OctTreeCell<X,E>::Y_MAX )
    {
      newmin.y( center.y() );
      newmax.y( oldmax.y() );
    }
    else
    {
      newmin.y( oldmin.y() );
      newmax.y( center.y() );
    }
    
    if( quadrant_flags & OctTreeCell<X,E>::Z_MAX )
    {
      newmin.z( center.z() );
      newmax.z( oldmax.z() );
    }
    else
    {
      newmin.z( oldmin.z() );
      newmax.z( center.z() );
    }
    
      // Put the box on the stack
    box_list[box_count++] = box;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Pop box/tree-node from stack and store in current_box
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
void OctTree<X,E>::pop_to_current()
{
  assert( box_count );
  box_count--;
  current_min = box_min[box_count];
  current_max = box_max[box_count];
  current_box = box_list[box_count];
  center = (current_min + current_max) * 0.5;
}

//-------------------------------------------------------------------------
// Purpose       : Allocate block of 8 oct-tree nodes for use as children
//                 of some (currently-leaf) node.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
OctTreeCell<X,E>* OctTree<X,E>::allocate_8()
{
    // Want to pop from end of page
  curr_page_end_ -= 8;
  
    // Any blocks available in current page?
  mem_pages_.last();
  if( curr_page_end_ < mem_pages_.get()  )
  {
      // Allocate new page
    mem_pages_.append( new OctTreeCell<X,E>[OCT_TREE_ALLOC_COUNT] );
    
      // Pop last block from new page
    mem_pages_.last();
    curr_page_end_ = mem_pages_.get() + OCT_TREE_ALLOC_COUNT - 8;
  }
  
  return curr_page_end_;
}


//-------------------------------------------------------------------------
// Purpose       : Recursively subdivide oct-tree nodes until
//                 one of three stop conditions are met:
//                   - min_nodes_ or less nodes in the box
//                   - size of child boxes will be smaller than
//                     min_box_size_ in all three dimensions
//                   - the nodes within this box are all within
//                     2*tolerance_ of each other.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/30/02
//-------------------------------------------------------------------------
template<class X, class E>
void OctTree<X,E>::split_node( OctTreeCell<X,E>* box, 
                             const CubitVector& min,
                             const CubitVector& max )
{
  assert( box->leaf() );
  CubitVector diag = max - min;
  diag *= 0.5; // diagonal size of child boxes
  if( (box->node_count() <= min_nodes_) ||  // must have more than min_nodes_
      (diag.x() < min_box_size_ &&          // child boxes will be smaller
       diag.y() < min_box_size_ &&          //   than min_box_size_ in all
       diag.z() < min_box_size_ ) )         //   three dimensions.
    return;
  
    // Check if all nodes are currently within 2*tolerance_
    // of each other.
  double tol2 = 2 * tolerance_;
  CubitVector node_min, node_max, junk;
  box->node_locations( node_min, node_max, junk );
  diag = node_max - node_min;
  if( diag.x() < tol2 &&
      diag.y() < tol2 &&
      diag.z() < tol2 )
    return;
  
    // Split the box
  CubitVector mid = (min + max) * 0.5;
  if( !box->split( mid, allocate_8() ) )
    return;
  
    // Recursively call split_node on all 8 new child boxes
  split_node( 
    box->child(OctTreeCell<X,E>::X_MIN|OctTreeCell<X,E>::Y_MIN|OctTreeCell<X,E>::Z_MIN),
    CubitVector( min.x(), min.y(), min.z() ),
    CubitVector( mid.x(), mid.y(), mid.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MIN|OctTreeCell<X,E>::Y_MIN|OctTreeCell<X,E>::Z_MAX),
    CubitVector( min.x(), min.y(), mid.z() ),
    CubitVector( mid.x(), mid.y(), max.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MIN|OctTreeCell<X,E>::Y_MAX|OctTreeCell<X,E>::Z_MIN),
    CubitVector( min.x(), mid.y(), min.z() ),
    CubitVector( mid.x(), max.y(), mid.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MIN|OctTreeCell<X,E>::Y_MAX|OctTreeCell<X,E>::Z_MAX),
    CubitVector( min.x(), mid.y(), mid.z() ),
    CubitVector( mid.x(), max.y(), max.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MAX|OctTreeCell<X,E>::Y_MIN|OctTreeCell<X,E>::Z_MIN),
    CubitVector( mid.x(), min.y(), min.z() ),
    CubitVector( max.x(), mid.y(), mid.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MAX|OctTreeCell<X,E>::Y_MIN|OctTreeCell<X,E>::Z_MAX),
    CubitVector( mid.x(), min.y(), mid.z() ),
    CubitVector( max.x(), mid.y(), max.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MAX|OctTreeCell<X,E>::Y_MAX|OctTreeCell<X,E>::Z_MIN),
    CubitVector( mid.x(), mid.y(), min.z() ),
    CubitVector( max.x(), max.y(), mid.z() ) );
  split_node( 
    box->child(OctTreeCell<X,E>::X_MAX|OctTreeCell<X,E>::Y_MAX|OctTreeCell<X,E>::Z_MAX),
    CubitVector( mid.x(), mid.y(), mid.z() ),
    CubitVector( max.x(), max.y(), max.z() ) );
}
  
  
template<class X, class E>
void OctTree<X,E>::tree_size( DLIList<int>& boxes, DLIList<int>& leaves )
{
  boxes.clean_out();
  boxes.reset();
  
  leaves.clean_out();
  leaves.reset();
  
  tree_size( root_, boxes, leaves );
}

template<class X, class E>
void OctTree<X,E>::tree_size( OctTreeCell<X,E>* box, DLIList<int>& boxes, DLIList<int>& leaves )
{
  if( boxes.is_at_end() )
  {
    boxes.append(1);
    boxes.step();
  }
  else
  {
    boxes.step();
    boxes.change_to( boxes.get() + 1 );
  }
  
  if( box->leaf() )
  {
    if( leaves.is_at_end() )
    {
      leaves.append(1);
      leaves.step();
    }
    else
    {
      leaves.step();
      leaves.change_to( leaves.get() + 1 );
    }
  }
  else
  {
    if( leaves.is_at_end() )
      leaves.append(0);
    leaves.step();

    for( int i = 0; i < 8; i++ )
      tree_size( box->child(i), boxes, leaves );
  }
  
  boxes.back();
  leaves.back();
}

#endif
