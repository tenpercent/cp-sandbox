///--------------------------------------------------------------------------
// Class Name:  RStarTreeNode
// Description: Node of R*Tree.
// The algorithm was taken from the following paper:
//	      Norbert Beckmann, H. Kriegel, R. Schnieder, and B. Seegar,
//              "The R*-tree: An Efficient and Robust Access Method
//              for Points and Rectangles", Proceedings of ACM SIGMOD
//              Int'l. Conf. on Management of Data, pp. 322-331, 1990.
// Creation Date: 7/21/02
// Owner:  David R. White
///--------------------------------------------------------------------------
#ifndef RSTARTREENODE_HPP
#define RSTARTREENODE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitBox.hpp"

template <class X> class DLIList;
const int UNSET_RSTARNODE = -1;
const int DATA_RSTARNODE = 0;
const int LEAF_RSTARNODE = 1;

template <class Y> class RStarTreeNode 
{
private:
    //Data.
  IttyBit markedFlag : 1;
    ///
    ///Generic mark flag.
    ///
  IttyBit distIsBox : 1;
  int myId;
  
    ///
    ///Mark if the distance measured is to the
    ///bounding box or the real data.  This is
    ///for the nearest neighbor search.
    ///
  RStarTreeNode<Y>** myChildrenNodes;
  int nextChildIndex;
  CubitBox *myBoundingBox;
  Y myData;
  int maxChildren, minChildren; //max/min number of children.
  RStarTreeNode<Y> *myParent;
  int myLevel;
    ///
    ///Level of node
    ///Level 0 equals data level.
    ///Level 1 equals leaf node level.
    ///Higher levels are non-leaf node levels.
    ///
  double myDist;
    ///
    ///used for nearest neigbhor searches...
    ///
  
  

    //Functions.
  RStarTreeNode<Y>* choose_sub_tree( RStarTreeNode<Y> *n, RStarTreeNode<Y> *e );
    ///
    /// Select a leaf node in which to place
    /// a new index entry e.  Recursive search the subtrees of n
    /// until n is a leaf node.
    ///

  CubitStatus overflow_treatment( RStarTreeNode<Y>* l,
                                  RStarTreeNode<Y>* e,
                                  RStarTreeNode<Y> *&ll,
                                  RStarTreeNode<Y> *root,
                                  RStarTreeNode<Y> *&new_root,
                                  int *overflow_flags, int levels);
    ///
    /// Determines if the node should be split or if nodes should be reinserted.
    /// If this is not the first time that the nodes have been overflown on this
    /// level (as determed by the overflow flags), then remove p items and reinsert
    /// them hoping to get a better box distribution.  Otherwise, split the node.
    /// We spend more time reinserting and less time splitting on average
    /// so the result is a better tree at basically the same cost.
    ///

  CubitStatus reinsert(RStarTreeNode<Y>* l,
                       RStarTreeNode<Y>* e,
                       RStarTreeNode<Y> *root,
                       RStarTreeNode<Y> *&new_root,
                       int *overflow_flags, int levels);
    ///
    /// Sorts the elements based on their distance to the centroid of the bounding rectangle.
    /// reinserts the top p nodes back into the tree (after removing them first of course.)
    ///

  
  CubitStatus adjust_tree(RStarTreeNode<Y> *l, RStarTreeNode<Y> *ll,
                          RStarTreeNode<Y> *root_node,
                          RStarTreeNode<Y> *&new_root,
                          int *overflow_flags,
                          int levels);
    ///
    /// Ascend from a leaf node L to the root, adjusting covering
    /// bounding boxes and propagating nodes splits as necesary.
    ///

  CubitStatus split_node( RStarTreeNode<Y> *l,
                          RStarTreeNode<Y> *e,
                          RStarTreeNode<Y> *&ll );
    ///
    /// Since element e won't fit into l, split l into two groups, l and ll.
    /// For the RTree, this is where most of the variations
    /// on implementations have occured.  Hopefully this is the best
    /// method for splitting.  Though I have seen some optimizations or variations
    /// on this that have claimed slightly better performance for given algorithms.
    ///
  double calc_overlap( RStarTreeNode<Y> *current,
                       RStarTreeNode<Y> *add_to);
    ///
    /// Calculates the current overlap between current's children nodes
    /// and add's to it the overlap that would occur from adding add_to to
    /// that list.
    ///
  
  double calc_enlargement( RStarTreeNode<Y> *current, RStarTreeNode<Y> *add_to );
  double calc_enlargement( CubitBox &current,
                           CubitBox &add_to );
    ///
    /// Calculate the enlargement required for increasing
    /// the bounding box of current so that it would encapsulate
    /// the bounding box of add_to.
    ///

#ifdef BOYD15
  CubitStatus pick_seeds(RStarTreeNode<Y> **input_list,
                         const int input_list_size,
                         RStarTreeNode<Y> *&seed_1,
                         RStarTreeNode<Y> *&seed_2);
    ///
    /// Picks the two starting children of the input list that
    /// are best for creating the two new groups of quadratic_split_node.
    ///
#endif

#ifdef BOYD15
  CubitStatus pick_next(DLIList <RStarTreeNode<Y>*> &remaining_nodes,
                        RStarTreeNode<Y>* group_1,
                        RStarTreeNode<Y>* group_2,
                        RStarTreeNode<Y>*& next,
                        CubitBoolean &add_to_group_1);
    ///
    /// picks one remainng entry for classification in a group.  Also
    /// indicates which group that should be by the add_to_group_1 flag.
    ///
#endif
  
  CubitStatus find_leaf( Y e,
                         CubitBox &e_box,
                         RStarTreeNode<Y> *t,
                         RStarTreeNode<Y> *&l );
    ///
    /// find the entry e in the tree t.  l
    /// is the entry that contains e. (e is a data member of l).
    /// The parent of l, is the leaf node containing it.
    ///

  CubitStatus condense_tree(RStarTreeNode<Y> *l,
                            RStarTreeNode<Y> *root,
                            RStarTreeNode<Y> *&new_root );
    ///
    ///  Given a leaf node l from which an entry has been deleted, eliminate
    ///  the node if it has too few entries and relocate its entries.  Propagate
    ///  node elimination upaward as necessary.  Adjust all covering rectangles
    ///  on the path to the root, making them smaller if possible.
    ///
    
  
  double volume(CubitBox &box);
    ///
    /// Calculate the volume of the box.
    ///
  double volume(RStarTreeNode<Y>* curr);
    ///
    /// Calculate the volume of the RStarTreeNode.
    ///

  static int sort_high_x(RStarTreeNode<Y> *&n_1,
                         RStarTreeNode<Y> *&n_2 );
  static int sort_high_y(RStarTreeNode<Y> *&n_1,
                         RStarTreeNode<Y> *&n_2 );
  static int sort_high_z(RStarTreeNode<Y> *&n_1,
                         RStarTreeNode<Y> *&n_2 );
  static int sort_low_x(RStarTreeNode<Y> *&n_1,
                        RStarTreeNode<Y> *&n_2);
  static int sort_low_y(RStarTreeNode<Y> *&n_1,
                        RStarTreeNode<Y> *&n_2);
  static int sort_low_z(RStarTreeNode<Y> *&n_1,
                        RStarTreeNode<Y> *&n_2);
    ///
    /// Functions to help sort the DLIList (thats why they are static).
    /// Returns -1 if n_1 is < than n_2 in that dimension, 0 if they
    /// are equivalent or 1.
    ///

  static int sort_center_distance( RStarTreeNode<Y> *&n_1,
                                   RStarTreeNode<Y> *&n_2 );
    ///
    /// Function for sorting during the reinsert function.
    /// Returns -1 if n_1->get_dist() is > than n_2, 0 if equal,
    /// else 1.
    ///

  

  double margin(CubitBox &bounding_box);
    ///
    /// Calculates the margin of the mounding box
    /// margin = ([xmax-xmin]+[ymax-ymin]+[zmax-zmin])*2^(d-1))
    ///

  CubitBox super_box(DLIList<RStarTreeNode<Y>*> &node_list);
    ///
    /// Finds the bounding box of all the nodes in the list.
    ///

public:
  RStarTreeNode(Y data, double tol, int max_children, int min_children );
    ///
    /// Constructor for data.
    ///
  RStarTreeNode(CubitBox &bounding_box, int max_children, int min_children);
    ///
    /// Constructor for typical parent node.
    ///
  
  ~RStarTreeNode();
    /// Constructor/Destructor

  void validate_tree(int print);
    /// Makes sure the children point to parents...

  CubitStatus insert(RStarTreeNode<Y> *e, RStarTreeNode<Y> *&new_root,
                     int *overflow_flags, int levels);
    ///
    /// Inserts the node e into the proper position of this.  Assumeing
    /// this is the root.  Sometimes the root will need to be split,
    /// so a new root may be returned.  If the new_root is not null,
    /// then the new_root must take the place of this...
    /// The overflow_flow flags variable needs to be initialized to
    /// be an aray of size equal to the number of levels + 1, and initialized
    /// to have values of zero.  The internal code uses this as to flag
    /// the levels as they get treated for overflow.
    ///

  void set_leaf_level(int r_type)
  {myLevel = r_type;}
    ///
    /// Set the type of level of the node.
    ///
  int get_leaf_level()
    {return myLevel;}
    ///
    /// get the level of the node.
    ///

  CubitBoolean is_leaf()
    {return (myLevel == LEAF_RSTARNODE)? CUBIT_TRUE : CUBIT_FALSE;}
    ///
    /// Determine if the RStarTreeNode is a leaf node.
    ///

  CubitBoolean is_data()
    {return (myLevel == DATA_RSTARNODE)? CUBIT_TRUE : CUBIT_FALSE;}
    ///
    /// Determine if the RStarTreeNode is a leaf node.
    ///

  Y get_data()
    {return myData;}
    ///
    /// Determine if the RStarTreeNode is a leaf node.
    ///
  
  void add_child(RStarTreeNode<Y>* child_node, CubitBoolean recalc_b_box);
    ///
    /// Add the child to the myChildrenNodes' list. Adds
    /// it to the next availabel spot.  Won't add if overflow will
    /// occur.
    ///
  
  CubitBoolean can_add();
    ///
    /// Tests if there is any space in the myChildrenNodes' list.
    ///

  int space_left();
    ///
    /// Returns the number of positions left in the myChildrenNode's list.
    ///

  int num_children()
    {return nextChildIndex;}
    ///
    /// Returns the number of children in the myChildrenNode's array.
    ///

  RStarTreeNode<Y>* get_child(int i)
    {return ((i < nextChildIndex) ? myChildrenNodes[i] : (RStarTreeNode<Y>*)NULL) ;}
  
    
  void flush(CubitBox &new_box);
    ///
    /// Clears out the myChildrenNodes by setting the array values
    /// to null.  Resets the counters and sets the range as the new_box.
    ///
  
#ifdef BOYD15
  void update_box( CubitBox &new_box );
    ///
    /// updates the nodes box with the new box.
    ///
#endif

  void recalc_b_box();
    ///
    /// recalculates the bounding box for the node. (won't do it if
    /// this is a data node...
    ///
  
  RStarTreeNode<Y> *get_parent()
    {return myParent;}
  void set_parent(RStarTreeNode<Y> *parent)
    {myParent = parent;}

  CubitBoolean remove_child(RStarTreeNode<Y> *child);
    ///
    /// Removes the child from the myChildrenNodes array and condenses
    /// the array.  decrements the number of children or increments the
    /// num positions available.
    ///

  CubitBoolean remove(Y e, RStarTreeNode<Y> *&new_root, CubitBoolean &delete_root);
    ///
    /// Removes the data e from the tree and cleans up the tree.  A new
    /// root node could be found from this.  Assume the root node is
    /// calling this function...
    ///

  CubitBox& bounding_box()
    {return *myBoundingBox;}
    ///
    /// Returns the bounding box of this node.
    ///

  double get_dist();
  void set_dist(double dis);
    ///
    /// This is used in several places to store some distance
    /// on the node.  It is used in sorting and overflow treatment.
    /// Basically at times the nodes need to be sorted according to
    /// some distance measurement.
    ///

  int dist_is_box();
  void set_dist_is_box(int val);
    ///
    /// Returns/sets the distIsBox flag used to determine if the
    /// stored distance is the distance between some object and this's
    /// bounding box or between some object and this's actual data.
    ///
  
      
};
template <class Y> inline int RStarTreeNode<Y>::dist_is_box()
{
  return distIsBox;
}
template <class Y> inline void RStarTreeNode<Y>::set_dist_is_box(int val)
{
  if ( val )
    distIsBox = 1;
  else
    distIsBox = 0;
}

template <class Y> inline double RStarTreeNode<Y>::get_dist()
{
  return myDist;
}
template <class Y> inline void RStarTreeNode<Y>::set_dist(double dist)
{
  myDist = dist;
}

//Calculate the volume of the box.
template <class Y> inline double RStarTreeNode<Y>::volume(CubitBox &box)
{
  return box.x_range()*box.y_range()*box.z_range();
}
//Calculate the volume of the RStarTreeNode.
template <class Y> inline double RStarTreeNode<Y>::volume(RStarTreeNode<Y>* curr)
{
  CubitBox box = curr->bounding_box();
  return box.x_range()*box.y_range()*box.z_range();
}
#if defined(TEMPLATE_DEFS_INCLUDED)
  #include "RStarTreeNode.cpp"
#endif

#endif

