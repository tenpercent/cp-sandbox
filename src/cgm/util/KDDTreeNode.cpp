//-----------------------------------------------------------------
//- Class:   KDDTreeNode
//- Author:  Kevin Albrecht 
//- Created: 13 May 2003
//- Updated: 8 Feb 2004
//-
//- Description:
//-   Node for dynamic version of the k-d tree, where k=3.
//-
//- References:
//-
//-   Hanan Samet.  Design and Analysis of Spatial Data Structures.
//-      Addison-Wesley, Reading, MA, 1990.
//-
//-   Jon Louis Bentley. Multidimensional binary search trees used
//-     for associative searching. In Communications of the ACM,
//-     18(9), pages 509-517, September 1975.
//-----------------------------------------------------------------

#if !defined(TEMPLATE_DEFS_INCLUDED) || defined(INCLUDED_FROM_KDD_TREE_NODE_HEADER)

//---------------------------------
// Include Files
//---------------------------------

#include "KDDTreeNode.hpp"
#include "DLIList.hpp"

//---------------------------
// Initialize Static Members
//---------------------------

#ifdef INLINE_TEMPLATES
  #define MY_INLINE inline
#else
  #define MY_INLINE
#endif

//- Constructor
template <class Y> MY_INLINE KDDTreeNode<Y>::KDDTreeNode
  ( Y aData, DIMENSION aDisc )
{
  parent = left = right = NULL;
  data = aData;

  boundingBox = data->bounding_box();
  x = boundingBox.center().x();
  y = boundingBox.center().y();
  z = boundingBox.center().z();

  myDist = CUBIT_DBL_MAX;
  myDistData = DD_SAFETY;
   
  myDisc = aDisc;

  valid = CUBIT_TRUE;
}

//- Destructor
template <class Y> MY_INLINE KDDTreeNode<Y>::~KDDTreeNode ()
{  
}

template <class Y> MY_INLINE KDDTreeNode<Y> *KDDTreeNode<Y>::get_child (DIRECTION dir) const
{
  if (dir == DIR_LEFT) return left;
  return right;
}

template <class Y> MY_INLINE void KDDTreeNode<Y>::set_child
  ( KDDTreeNode<Y> *node,
    DIRECTION dir )
{
  if (dir == DIR_LEFT)
  {
    left = node;
  }
  else
  {
    right = node;
  }
}

template <class Y> MY_INLINE DIMENSION KDDTreeNode<Y>::next_disc () const
{
  switch (myDisc)
  {
    case DIMX: return DIMY;
    case DIMY: return DIMZ;
    default:   return DIMX;
  }
}

//- The KD_COMPARE function as defined by Samet
template <class Y> MY_INLINE DIRECTION KDDTreeNode<Y>::compare (KDDTreeNode<Y> *Q) const
{
  if (Q->myDisc == DIMX)
  {
    if (x <= Q->x) return DIR_LEFT;
    else return DIR_RIGHT;
  }
  else if (Q->myDisc == DIMY)
  {
    if (y <= Q->y) return DIR_LEFT;
    else return DIR_RIGHT;
  }
  else if (z <= Q->z)
  {
    return DIR_LEFT;
  }
  else
  {
    return DIR_RIGHT;
  }
}

//- The KD_COMPARE function as defined by Samet
template <class Y> MY_INLINE DIRECTION KDDTreeNode<Y>::compare_with_equality (KDDTreeNode<Y> *Q) const
{
  if (Q->myDisc == DIMX)
  {
    if (x < Q->x) return DIR_LEFT;
    else if (x == Q->x) return DIR_EITHER;
    else return DIR_RIGHT;
  }
  else if (Q->myDisc == DIMY)
  {
    if (y < Q->y) return DIR_LEFT;
    else if (y == Q->y) return DIR_EITHER;
    else return DIR_RIGHT;
  }
  else if (z < Q->z)
  {
    return DIR_LEFT;
  }
  else if (z == Q->z)
  {
    return DIR_EITHER;
  }
  else
  {
    return DIR_RIGHT;
  }
}

#endif
