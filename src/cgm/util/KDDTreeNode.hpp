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

#ifndef KDDTREENODE_HPP
#define KDDTREENODE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitBox.hpp"

//- constants that show the direction of discrimination
enum DIMENSION
{
  DIMX, DIMY, DIMZ
};

//- constants used to determine left vs. right
enum DIRECTION
{
  DIR_NULL, DIR_EITHER, DIR_LEFT, DIR_RIGHT
};

enum DISTDATA
{
  DD_SAFETY, DD_LEAF
};

template <class X> class DLIList;

template <class Y> class KDDTreeNode 
{
  private:
    DIMENSION myDisc;         // dimension in which this node discriminates
    double myDist;            // used for nearest neigbhor searches
    DISTDATA myDistData;      // also used for nearest neigbhor searches

  public:
    KDDTreeNode<Y> *parent; // parent of this node
    KDDTreeNode<Y> *left;   // left subtree
    KDDTreeNode<Y> *right;  // right subtree
    double x, y, z;         // coordinates of center point
    Y data;                 // data with bounding box centered on this node
    CubitBoolean valid;     // false if this node should be deleted

    CubitBox safetyBox;   // bounding box of myData and the safetyBox's of children
    CubitBox boundingBox; // bounding box of myData

    //- Constructor/Destructor
    KDDTreeNode (Y aData, DIMENSION aDisc = DIMX);
    ~KDDTreeNode ();

    //- Get/set the child in the direction of "dir"
    KDDTreeNode<Y> *get_child (DIRECTION dir) const;
    void set_child (KDDTreeNode<Y> *node, DIRECTION dir);

    //- Get/set the "myDist" value
    double get_dist () const { return myDist; }
    void set_dist (double dist) { myDist = dist; }

    //- Get/set the "myDistData" value
    DISTDATA get_dist_data () const { return myDistData; }
    void set_dist_data (DISTDATA distData) { myDistData = distData; }

    //- Get/set the discriminating dimension of this node
    DIMENSION get_disc () const { return myDisc; }
    void set_disc (DIMENSION dim) { myDisc = dim; }

    //- Get the discriminator for the depth below this node
    DIMENSION next_disc () const;

#ifdef BOYD15
    //- Get the center point of this node
    CubitVector get_center () const;
#endif

    //- The KD_COMPARE function as defined by Samet
    DIRECTION compare (KDDTreeNode<Y> *Q) const;
    DIRECTION compare_with_equality (KDDTreeNode<Y> *Q) const;
};

#ifdef TEMPLATE_DEFS_INCLUDED
#  define INCLUDED_FROM_KDD_TREE_NODE_HEADER
#  include "KDDTreeNode.cpp"
#  undef INCLUDED_FROM_KDD_TREE_NODE_HEADER
#endif

#endif
