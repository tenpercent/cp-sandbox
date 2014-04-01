//-----------------------------------------------------------------
//- Class:   KDDTree
//- Author:  Kevin Albrecht 
//- Created: 13 May 2003
//- Updated: 8 Feb 2004
//-
//- Description:
//-   Dynamic version of the k-d tree, where k=3.
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

#ifndef KDDTREE_HPP
#define KDDTREE_HPP

// this value is a percentage
#define KDDTREE_DEFAULT_SELF_BALANCING_DELETION_TOLERANCE 0.0

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "KDDTreeNode.hpp"
#include "PriorityQueue.hpp"
#include "AbstractTree.hpp"

class CubitBox;
class CubitVector;

template <class X> class DLIList;
template <class Y> class KDDTreeNode;

template <class Z> class KDDTree : public AbstractTree<Z>
{
  private:
    CubitBoolean mySelfBalancingOn;
    double myTolerance;
    DLIList<KDDTreeNode<Z>*> myAddList;  // nodes added but not on tree
    KDDTreeNode<Z> *myDeepestLeaf;       // used by "recursiveRemove" to hold the deepest leaf node
    int myMarkedNodes;                   // number of nodes marked for deletion
    double myDeletionTolerance;          // the deletion tolerance (floating point percentage between 0 and 1)
    CubitBoolean myRandomOn;             // randomness switch

    double min_dist_sq (CubitVector &q, CubitBox &b_box);

    //- Find the depth of the tree recursively
    void recursive_find_max_height (KDDTreeNode<Z> *the_root, int depth, int *maxdepth);

    //- Balance the tree recursively
    KDDTreeNode<Z> *recursive_balance (DIMENSION dim, int left, int right,
                                      KDDTreeNode<Z>* array[], KDDTreeNode<Z>* parent);

    //- Recursively find members intersecting this range box (called by "find")
    void recursive_find ( KDDTreeNode<Z> *rect_tree, const CubitBox &range_box,
                          DLIList <Z> &range_members );

    //- Return a pointer to the node containing the specified data
    KDDTreeNode<Z> *find_node_containing_data (KDDTreeNode<Z> *subtreeRoot, Z data);

    //- Find the node with the minimum value in the D dimension (used by
    //-  "recursiveRemove")
    KDDTreeNode<Z> *find_minimum_node_of_dimension (KDDTreeNode<Z> *P, DIMENSION D);

    //- Insert the node on the tree
    CubitStatus insert_node (KDDTreeNode<Z>* P);

    //- Immediately put all nodes on list onto the tree
    CubitStatus dump_list ();

    //- Rearrange the array around the median point
    int modifind (DIMENSION dim, int left, int right, KDDTreeNode<Z>* array[]);

  public:
    DLIList<KDDTreeNode<Z>*> myNodeList; // nodes on tree
    KDDTreeNode<Z> *root;

    //- Constructors/Destructor
    KDDTree (double tol = GEOMETRY_RESABS, CubitBoolean selfBalancingOn = CUBIT_TRUE,
             double selfBalancingDeletionTolerance = KDDTREE_DEFAULT_SELF_BALANCING_DELETION_TOLERANCE,
             CubitBoolean randomOn = CUBIT_FALSE);
    ~KDDTree();

    //- Methods for finding the maximum height of the tree and the size
    //- of the tree
    int find_max_height ();

    //- Sets/gets the tolerance used for the bounding box overlap test,
    //- which is used during the range search
    void set_tol (double tol) { myTolerance = tol; }
    double get_tol () { return myTolerance; }

    //- Add a node with the data to the list
    CubitStatus add (Z data);

    //- Remove the data member's entry in the tree. Returns CUBIT_TRUE
    //- if item removed, CUBIT_FALSE if item not in tree
    CubitBoolean remove (Z data);

    //- Balance the tree manually
    //- Note: this is also called by the find method when self-balancing is on
    CubitStatus balance ();

    //- Find members intersecting this range box
    CubitStatus find (const CubitBox &range_box, DLIList <Z> &range_members);

    //- Distance method used by k_nearest_neighbor
    typedef double (*DistSqFunc)(CubitVector &a, Z& b);

    //- Find the k nearest neighbors
    CubitStatus k_nearest_neighbor (CubitVector &q, int k, double &closest_dist,
                                    DLIList<Z> &nearest_neighbors,
                                    DistSqFunc dist_sq_point_data
                                   );

    //- Function to use with the priority queue for sorting.
    static bool less_than_func (KDDTreeNode<Z> *&node_a,
                                KDDTreeNode<Z> *&node_b);
};

#ifdef TEMPLATE_DEFS_INCLUDED
#  define INCLUDED_FROM_KDD_TREE_HEADER
#  include "KDDTree.cpp"
#  undef INCLUDED_FROM_KDD_TREE_HEADER
#endif

#endif
