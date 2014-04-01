//-------------------------------------------------------------------------
// Filename      : OctTreeCell.cpp
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

#if !defined(TEMPLATE_DEFS_INCLUDED) || defined(INCLUDED_FROM_OCT_TREE_CELL_HPP)

#include "OctTreeCell.hpp"
#include "CubitVector.hpp"

template <class X, class E>
OctTreeCell<X,E>::OctTreeCell( OctTreeEntry<X,E> array[], int size )
  : node_count_(size)
{
  head_or_children = 0;
  OctTreeEntry<X,E>** prev_ptr_ptr = reinterpret_cast<OctTreeEntry<X,E>**>(&head_or_children);
  for( int i = 0; i < node_count_; i++ )
  {
    OctTreeEntry<X,E>* entry = array + i;
    *prev_ptr_ptr = entry;
    prev_ptr_ptr = &(entry->next);
  }
  *prev_ptr_ptr = 0;
}

template <class X, class E>
void OctTreeCell<X,E>::append_nodes( DLIList<X*>& list )
{
  OctTreeEntry<X,E>* entry = node_count_ ? reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children) : 0;
  for( ; entry; entry = entry->next )
    list.append( entry->node );
}
  
template <class X, class E>
void OctTreeCell<X,E>::append_all_nodes( DLIList<X*>& list )
{
  if( node_count_ )
  {
    append_nodes( list );
  }
  else if( reinterpret_cast<OctTreeCell<X,E>*>(head_or_children) )
  {
    for( int i = 0; i < 8; i++ )
      reinterpret_cast<OctTreeCell<X,E>*>(head_or_children)[i].append_all_nodes( list );
  }
}

template <class X, class E>
void OctTreeCell<X,E>::add_nodes( OctTreeEntry<X,E>* node )
{
  assert( leaf() );
  
  OctTreeEntry<X,E>* last = node;
  int count = 1;
  while( last->next ) 
  {
    count++;
    last = last->next;
  }
  
  last->next = reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children);
  head_or_children = reinterpret_cast<void*>(node);
  node_count_ += count;
}

template <class X, class E>
OctTreeCell<X,E>* OctTreeCell<X,E>::child( int quadrant )
{
  assert( quadrant >= 0 && quadrant < 8 && !node_count_ );
  if( reinterpret_cast<OctTreeCell<X,E>*>(head_or_children) )
    return reinterpret_cast<OctTreeCell<X,E>*>(head_or_children) + quadrant;
  else
    return 0;
}

template <class X, class E>
bool OctTreeCell<X,E>::split( const CubitVector& my_center,
                            OctTreeCell<X,E>* storage )
{
  assert( leaf() );
  bool result = false;
  
  if( node_count_ > 3 ) 
  {
      // Data is a union.  Make sure you never change
      // the order of the following two lines!
    OctTreeEntry<X,E>* node = reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children);
    head_or_children = reinterpret_cast<void*>(storage);

    while( node )
    {
      CubitVector coords = E::coordinates(node->node);
      int x = coords.x() < my_center.x() ? X_MIN : X_MAX;
      int y = coords.y() < my_center.y() ? Y_MIN : Y_MAX;
      int z = coords.z() < my_center.z() ? Z_MIN : Z_MAX;
      OctTreeEntry<X,E>* next = node->next;
      node->next = 0;
      reinterpret_cast<OctTreeCell<X,E>*>(head_or_children)[x|y|z].add_nodes( node );
      node = next;
    }
    result = true;
    node_count_ = 0;
  }
  
  return result;
}  
  
template <class X, class E>
void OctTreeCell<X,E>::node_locations( CubitVector& box_min,
                                     CubitVector& box_max,
                                     CubitVector& average,
                                     CubitVector* std_dev )
{
  if( !node_count_ )
    return;


  double x_min, y_min, z_min;
  double x_max, y_max, z_max;
  double x_avg, y_avg, z_avg;
  CubitVector coords = E::coordinates( reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children)->node );
  x_min = x_max = x_avg = coords.x();
  y_min = y_max = y_avg = coords.y();
  z_min = z_max = z_avg = coords.z();
  
  for( OctTreeEntry<X,E>* entry = reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children)->next; entry; entry = entry->next )
  {
    X* node = entry->node;
    coords = E::coordinates( node );

    double x = coords.x();
    if( x < x_min ) x_min = x;
    else if ( x > x_max ) x_max = x;
    x_avg += x;
    
    double y = coords.y();
    if( y < y_min ) y_min = y;
    else if ( y > y_max ) y_max = y;
    y_avg += y;
    
    double z = coords.z();
    if( z < z_min ) z_min = z;
    else if ( z > z_max ) z_max = z;
    z_avg += z;
  }
  
  box_min.set( x_min, y_min, z_min );
  box_max.set( x_max, y_max, z_max );
  
  double inv = 1.0 / node_count_;
  x_avg *= inv;
  y_avg *= inv;
  z_avg *= inv;
  average.set( x_avg, y_avg, z_avg );
  
  if( std_dev )
  {
    X* hn = reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children)->node;
    coords = E::coordinates(hn);
    double x_dev = (coords.x() - x_avg) * (coords.x() - x_avg);
    double y_dev = (coords.y() - y_avg) * (coords.y() - y_avg);
    double z_dev = (coords.z() - z_avg) * (coords.z() - z_avg);

    for( OctTreeEntry<X,E>* entry = reinterpret_cast<OctTreeEntry<X,E>*>(head_or_children)->next; entry; entry = entry->next )
    {
      X* node = entry->node;
      coords = E::coordinates(node);

      double x = coords.x();
      double dx = x - x_avg;
      x_dev += dx * dx;

      double y = coords.y();
      double dy = y - y_avg;
      y_dev += dy * dy;

      double z = coords.z();
      double dz = z - z_avg;
      z_dev += dz * dz;
    }

    inv = 1.0 / ( node_count_ - 1 );
    x_dev = sqrt(x_dev * inv);
    y_dev = sqrt(y_dev * inv);
    z_dev = sqrt(z_dev * inv);

    std_dev->set( x_dev, y_dev, z_dev );
  }
}


#endif
