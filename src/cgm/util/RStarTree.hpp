//---------------------------------------------------------------------------
// Class Name:  RStarTree
// Description: Rectangle tree.  A tree for multidimensional access methods.
// The algorithm was taken from the following paper:
//	      Norbert Beckmann, H. Kriegel, R. Schnieder, and B. Seegar,
//              "The R*-tree: An Efficient and Robust Access Method
//              for Points and Rectangles", Proceedings of ACM SIGMOD
//              Int'l. Conf. on Management of Data, pp. 322-331, 1990.
// Creation Date: 7/20/02
// Owner:  David R. White
//---------------------------------------------------------------------------
#ifndef RSTARTREE_HPP
#define RSTARTREE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"

class CubitBox;
class CubitVector;

template <class X> class DLIList;
template <class Y> class RStarTreeNode;

template <class Z> class RStarTree 
{
private:
  double myTolerance;
  int maxChildren;
  int minChildren;
  RStarTreeNode<Z> *myRoot;
  CubitStatus recursive_find(RStarTreeNode<Z> *rect_tree,
                             const CubitBox &range_box,
                             DLIList <Z> &range_members);
    //- recurses down the rtree to find all all the data members
    //- that fall in the range of the range_box.

  void to_list(DLIList <RStarTreeNode<Z>*> &member_list,
               RStarTreeNode<Z> *top);
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
  
  RStarTree(double tol = GEOMETRY_RESABS);
  RStarTree(double tol, int max_c, int min_c);
  ~RStarTree();
    //- Constructor/Destructor

#ifdef BOYD15
  CubitStatus add(Z data);
    //- Adds the data member to the RStarTree.
#endif

#ifdef BOYD15
  CubitStatus find( const CubitBox &range_box, DLIList <Z> &range_members);
    //- searches the range tree for members that intersect this range box
    //- within the tolerance.

  CubitBoolean remove(Z data );
    //- Remove the data member's entry in the rectangle tree.
    //- Returns CUBIT_TRUE if item removed.  FALSE if item not
    //- in tree.
#endif

  typedef double (*DistSqFunc)(CubitVector &a, Z& b);
#ifdef BOYD15
  CubitStatus k_nearest_neighbor(CubitVector &q,
                                 int k,
                                 double &closest_dist,
                                 DLIList<Z> &nearest_neighbors,
                                 DistSqFunc dist_sq_point_data);
#endif
  
  void set_tol(double tol)
    {myTolerance = tol;}
  double get_tol()
    {return myTolerance;}
    //- Sets/Gets the tolerance used for the bounding box overlap test,
    //- which is used during the range search.

#ifdef BOYD15
  static bool less_than_func(RStarTreeNode<Z> *&node_a,
                             RStarTreeNode<Z> *&node_b);
    //- Function to use with the priority queue for sorting.
#endif
  
};
#if defined(TEMPLATE_DEFS_INCLUDED)
 #include "RStarTree.cpp"
#endif

#endif


