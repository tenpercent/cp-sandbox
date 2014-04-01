//---------------------------------------------------------------------------
// Class Name:  RTreeNode
// Description: Node of Rectangle tree.  Contians many of the
//              required functions for building the tree and traversing it.
// The algorithm was taken from the following paper:
//	      Guttman, A., "R-Trees: A Dynamic Index Structure for
//              Spatial Searching", Proceedings of the SIGMOD
//              Conference, Boston, June 1984, p. 47-57.
// Creation Date: 3/13/02
// Owner:  David R. White
//---------------------------------------------------------------------------
#ifndef RTREENODE_HPP
#define RTREENODE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitBox.hpp"

template <class X> class DLIList;
const int UNSET_RNODE = -1;
const int DATA_RNODE = 0;
const int LEAF_RNODE = 1;

template <class Y> class RTreeNode 
{
private:
    //Data.
  IttyBit markedFlag : 1;
  IttyBit distIsBox : 1;
    //Generic mark flag.
  RTreeNode<Y>** myChildrenNodes;
  int nextChildIndex;
  CubitBox *myBoundingBox;
  Y myData;
  int maxChildren, minChildren; //max/min number of children.
  RTreeNode<Y> *myParent;
  int myLevel; //Level of node
    //Level 0 equals data level.
    //Level 1 equals leaf node level.
    //Higher levels are non-leaf node levels.
  double myDist; //used for nearest neigbhor searches...
  

    //Functions.
  RTreeNode<Y>* choose_leaf( RTreeNode<Y> *n, RTreeNode<Y> *e );
    //- Select a leaf node in which to place
    //- a new index entry e.  Recursive search the subtrees of n
    //- until n is a leaf node.

  CubitStatus adjust_tree(RTreeNode<Y> *l, RTreeNode<Y> *ll,
                          RTreeNode<Y> *root_node,
                          RTreeNode<Y> *&new_root);
    //- Ascend from a leaf node L to the root, adjusting covering
    //- bounding boxes and propagating nodes splits as necesary.

  CubitStatus quadratic_split_node( RTreeNode<Y> *l,
                                    RTreeNode<Y> *e,
                                    RTreeNode<Y> *&ll );
    //- Since element e won't fit into l, split l into two groups, l and ll.
    //- For the RTree, this is where most of the variations
    //- on implementations have occured.  The best method proposed by Guttman
    //- was a quadratic split, which I'll implement here.  The Rstar tree
    //- did some slightly more complicated things which I might try later.
  
  
  double calc_enlargement( RTreeNode<Y> *current, RTreeNode<Y> *add_to );
  double calc_enlargement( CubitBox &current,
                           CubitBox &add_to );
    //- Calculate the enlargement required for increasing
    //- the bounding box of current so that it would encapsulate
    //- the bounding box of add_to.

  CubitStatus pick_seeds(RTreeNode<Y> **input_list,
                         const int input_list_size,
                         RTreeNode<Y> *&seed_1,
                         RTreeNode<Y> *&seed_2);
    //- Picks the two starting children of the input list that
    //- are best for creating the two new groups of quadratic_split_node.

  CubitStatus pick_next(DLIList <RTreeNode<Y>*> &remaining_nodes,
                        RTreeNode<Y>* group_1,
                        RTreeNode<Y>* group_2,
                        RTreeNode<Y>*& next,
                        CubitBoolean &add_to_group_1);
    //- picks one remainng entry for classification in a group.  Also
    //- indicates which group that should be by the add_to_group_1 flag.
  
  CubitStatus find_leaf( Y e,
                         CubitBox &e_box,
                         RTreeNode<Y> *t,
                         RTreeNode<Y> *&l );
    //- find the entry e in the tree t.  l
    //- is the entry that contains e. (e is a data member of l).
    //- The parent of l, is the leaf node containing it.

  CubitStatus condense_tree(RTreeNode<Y> *l,
                            RTreeNode<Y> *root,
                            RTreeNode<Y> *&new_root );
  //-  Given a leaf node l from which an entry has been deleted, eliminate
  //-  the node if it has too few entries and relocate its entries.  Propagate
  //-  node elimination upaward as necessary.  Adjust all covering rectangles
  //-  on the path to the root, making them smaller if possible.
    
  
  double volume(CubitBox &box);
  //- Calculate the volume of the box.
  double volume(RTreeNode<Y>* curr);
  //- Calculate the volume of the RTreeNode.


public:
  RTreeNode(Y data, double tol, int max_children, int min_children );
  RTreeNode(CubitBox &bound_box, int max_children, int min_children);
  
  ~RTreeNode();
    //- Constructor/Destructor

  CubitStatus insert(RTreeNode<Y> *e, RTreeNode<Y> *&new_root);
    //- Inserts the node e into the proper position of this.  Assumeing
    //- this is the root.  Sometimes the root will need to be split,
    //- so a new root may be returned.  If the new_root is not null,
    //- then the new_root must take the place of this...

  void set_leaf_level(int r_type)
  {myLevel = r_type;}
    //- Set the type of level of the node.
  int get_leaf_level()
    {return myLevel;}
    //- get the level of the node.

  CubitBoolean is_leaf()
    {return (myLevel == LEAF_RNODE)? CUBIT_TRUE : CUBIT_FALSE;}
    //- Determine if the RTreeNode is a leaf node.

  CubitBoolean is_data()
    {return (myLevel == DATA_RNODE)? CUBIT_TRUE : CUBIT_FALSE;}
    //- Determine if the RTreeNode is a leaf node.

  Y get_data()
    {return myData;}
    //- Determine if the RTreeNode is a leaf node.
  
  void add_child(RTreeNode<Y>* child_node, CubitBoolean recalc_bound_box);
    //- Add the child to the myChildrenNodes' list. Adds
    //- it to the next availabel spot.  Won't add if overflow will
    //- occur.
  
  CubitBoolean can_add();
    //- Tests if there is any space in the myChildrenNodes' list.

  int space_left();
    //- Returns the number of positions left in the myChildrenNode's list.

  int num_children()
    {return nextChildIndex;}
    //- Returns the number of children in the myChildrenNode's array.

  RTreeNode<Y>* get_child(int i)
    {return ( (i < nextChildIndex) ? myChildrenNodes[i] : static_cast< RTreeNode<Y>* >(NULL) );}
  
    
  void flush(CubitBox &new_box);
    //- Clears out the myChildrenNodes by setting the array values
    //- to null.  Resets the counters and sets the range as the new_box.
  
#ifdef BOYD15
  void update_box( CubitBox &new_box );
    //- updates the nodes box with the new box.
#endif

  void recalc_b_box();
    //- recalculates the bounding box for the node. (won't do it if
    //- this is a data node...
  
  RTreeNode<Y> *get_parent()
    {return myParent;}
  void set_parent(RTreeNode<Y> *parent)
    {myParent = parent;}

  CubitBoolean remove_child(RTreeNode<Y> *child);
    //- removes the child from the myChildrenNodes array and condenses
    //- the array.  decrements the number of children or increments the
    //- num positions available.

  CubitBoolean remove(Y e, RTreeNode<Y> *&new_root, CubitBoolean &delete_root);
    //- removes the data e from the tree and cleans up the tree.  A new
    //- root node could be found from this.  Assume the root node is
    //- calling this function...

  CubitBox& bounding_box()
    {return *myBoundingBox;}
    //- Returns the bounding box of this node.

  double get_dist();
  void set_dist(double dis);

  int dist_is_box();
  void set_dist_is_box(int val);
  
      
};
template <class Y> inline int RTreeNode<Y>::dist_is_box()
{
  return distIsBox;
}
template <class Y> inline void RTreeNode<Y>::set_dist_is_box(int val)
{
  if ( val )
    distIsBox = 1;
  else
    distIsBox = 0;
}

template <class Y> inline double RTreeNode<Y>::get_dist()
{
  return myDist;
}
template <class Y> inline void RTreeNode<Y>::set_dist(double dist)
{
  myDist = dist;
}

//Calculate the volume of the box.
template <class Y> inline double RTreeNode<Y>::volume(CubitBox &box)
{
  return box.x_range()*box.y_range()*box.z_range();
}
//Calculate the volume of the RTreeNode.
template <class Y> inline double RTreeNode<Y>::volume(RTreeNode<Y>* curr)
{
  CubitBox box = curr->bounding_box();
  return box.x_range()*box.y_range()*box.z_range();
}
#if defined(TEMPLATE_DEFS_INCLUDED)
  #include "RTreeNode.cpp"
#endif

#endif

