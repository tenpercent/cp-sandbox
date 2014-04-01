//---------------------------------------------------------------------------
// Class Name:  RStarTree
// Description: Rectangle tree.  Multidimensional access method (efficient
//              method to find ranges of boxes.
// The algorithm was taken from the following paper:
//	      Norbert Beckmann, H. Kriegel, R. Schnieder, and B. Seegar,
//              "The R*-tree: An Efficient and Robust Access Method
//              for Points and Rectangles", Proceedings of ACM SIGMOD
//              Int'l. Conf. on Management of Data, pp. 322-331, 1990.
// Creation Date: 7/20/02
// Owner:  David R. White
//---------------------------------------------------------------------------

//---------------------------------
//Include Files
//---------------------------------
#include "RStarTree.hpp"
#include "RStarTreeNode.hpp"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "PriorityQueue.hpp"
//---------------------------
//Initialize Static Members
//---------------------------

#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

template <class Z> MY_INLINE RStarTree<Z>::~RStarTree()
{
  if ( myRoot != NULL )
  {
      //Go through and get all the children in a list.
    DLIList <RStarTreeNode<Z>*> to_delete;
    to_list(to_delete, myRoot);
    int ii;
    for(ii = to_delete.size(); ii > 0; ii-- )
      delete to_delete.pop();
    delete myRoot;
  }
}
template <class Z> MY_INLINE void RStarTree<Z>::to_list(DLIList <RStarTreeNode<Z>*> &member_list,
                                          RStarTreeNode<Z> *top)
{
    //Get the children of the top into the list.
  int ii;
  RStarTreeNode <Z> *curr_node;
  for ( ii = 0; ii < top->num_children(); ii++ )
  {
    curr_node = top->get_child(ii);
    member_list.append(curr_node);
      //don't go below the bottom level...
    if ( curr_node->get_leaf_level() == 0 )
      continue;
    to_list(member_list, curr_node);
  }
  return;
}

template <class Z> MY_INLINE CubitStatus RStarTree<Z>::recursive_find(RStarTreeNode<Z> *rect_tree,
                                                        const CubitBox &range_box,
                                                        DLIList <Z> &range_members )
{
  CubitBox rect_box = rect_tree->bounding_box();
  if ( !range_box.overlap(myTolerance, rect_box ) )
    return CUBIT_SUCCESS;

    //Now see if this is a data member.  If it is, append the data to the
    //list.
  if (rect_tree->is_data() )
  {
    range_members.append(rect_tree->get_data());
    return CUBIT_SUCCESS;
  }
    //Now if this is anything else we need to keep iterating...
  int loop_size = rect_tree->num_children();
    //We are doing a depth-first search of the tree.  Not
    //all branches will need to be followed since they won't
    //all overlap...
  int ii;
  RStarTreeNode<Z> *curr_node;
  CubitStatus stat;
  for ( ii = 0; ii < loop_size; ii++ )
  {
    curr_node = rect_tree->get_child(ii);
    if ( curr_node == NULL )
    {
      PRINT_ERROR("Problems finding boxes in range.\n");
      assert(curr_node != NULL);
      return CUBIT_FAILURE;
    }
    stat = recursive_find(curr_node, range_box, range_members);
    if ( stat != CUBIT_SUCCESS )
      return stat;
  }
  
  return CUBIT_SUCCESS;
}
//--------------------------------------------------------------------------
//Algorithm: min_dist_sq
//Description:  Finds the minimum distance squared between the given
//              point and the box. If the point is on or in the box, the
//              min distance is zero.
//--------------------------------------------------------------------------
template <class Z> MY_INLINE
double RStarTree<Z>::min_dist_sq(CubitVector &q,
                             CubitBox &b_box)
{
  CubitVector b_min, b_max;
  b_min = b_box.minimum();
  b_max = b_box.maximum();
  double dist;
  CubitVector r;

  if ( q.x() < b_min.x() )
    r.x(b_min.x());
  else if ( q.x() > b_max.x() )
    r.x(b_max.x());
  else
    r.x(q.x());
  
  if ( q.y() < b_min.y() )
    r.y(b_min.y());
  else if ( q.y() > b_max.y() )
    r.y(b_max.y());
  else
    r.y(q.y());
  
  if ( q.z() < b_min.z() )
    r.z(b_min.z());
  else if ( q.z() > b_max.z() )
    r.z(b_max.z());
  else
    r.z(q.z());
  
  dist = (q-r).length_squared();

  return dist;
}

    
    
