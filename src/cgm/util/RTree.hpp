//---------------------------------------------------------------------------
// Class Name:  RTree
// Description: Rectangle tree.  A tree for multidimensional access methods.
// The algorithm was taken from the following paper:
//	      Guttman, A., "R-Trees: A Dynamic Index Structure for
//              Spatial Searching", Proceedings of the SIGMOD
//              Conference, Boston, June 1984, p. 47-57.
// Creation Date: 3/29/02
// Owner:  David R. White
//---------------------------------------------------------------------------
#ifndef RTREE_HPP
#define RTREE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "AbstractTree.hpp"

class CubitBox;
class CubitVector;

template <class X> class DLIList;
template <class Y> class RTreeNode;

template <class Z> class RTree : public AbstractTree<Z>
{
private:
  double myTolerance;
  int maxChildren;
  int minChildren;
  RTreeNode<Z> *myRoot;
  CubitStatus recursive_find(RTreeNode<Z> *rect_tree,
                             const CubitBox &range_box,
                             DLIList <Z> &range_members);
    //- recurses down the rtree to find all all the data members
    //- that fall in the range of the range_box.

  void to_list(DLIList <RTreeNode<Z>*> &member_list,
               RTreeNode<Z> *top);
    //- converts the tree under top to a list.  Note, top is NOT in the list
    //- at the end.
  double min_dist_sq(CubitVector &q,
                     CubitBox &b_box);
    //- Calculates the distance between a point and a box, this is
//- for the nearest neighbor point search.

  
   class LessThan
   {
     public:
     bool operator()(const Z& x, const Z& y) const
       {
         return x->get_dist() < y->get_dist(); 
       }
     
   };


  
public:
  RTree(double tol = GEOMETRY_RESABS);
  RTree(double tol, int max_c, int min_c);
  ~RTree();
    //- Constructor/Destructor

  CubitStatus add(Z data);
    //- Adds the data member to the RTree.

  CubitStatus find( const CubitBox &range_box, DLIList <Z> &range_members);
    //- searches the range tree for members that intersect this range box
    //- within the tolerance.

  CubitBoolean remove(Z data );
    //- Remove the data member's entry in the rectangle tree.
    //- Returns CUBIT_TRUE if item removed.  FALSE if item not
    //- in tree.

  typedef double (*DistSqFunc)(CubitVector &a, Z& b);
  CubitStatus k_nearest_neighbor(CubitVector &q,
                                 int k,
                                 double &closest_dist,
                                 DLIList<Z> &nearest_neighbors,
                                 DistSqFunc dist_sq_point_data);
  
  void set_tol(double tol)
    {myTolerance = tol;}
  double get_tol()
    {return myTolerance;}
    //- Sets/Gets the tolerance used for the bounding box overlap test,
    //- which is used during the range search.

  static bool less_than_func(RTreeNode<Z> *&node_a,
                             RTreeNode<Z> *&node_b);
    //- Function to use with the priority queue for sorting.
  
};
#if defined(TEMPLATE_DEFS_INCLUDED)
 #include "RTree.cpp"
#endif

#endif

