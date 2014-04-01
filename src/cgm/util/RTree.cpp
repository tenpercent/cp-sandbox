//---------------------------------------------------------------------------
// Class Name:  RTree
// Description: Rectangle tree.  Multidimensional access method (efficient
//              method to find ranges of boxes.
// The algorithm was taken from the following paper:
//	      Guttman, A., "R-Trees: A Dynamic Index Structure for
//              Spatial Searching", Proceedings of the SIGMOD
//              Conference, Boston, June 1984, p. 47-57.
// Creation Date: 3/29/02
// Owner:  David R. White
//---------------------------------------------------------------------------

//---------------------------------
//Include Files
//---------------------------------
#include "RTree.hpp"
#include "RTreeNode.hpp"
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
template <class Z> MY_INLINE RTree<Z>::RTree (double tol)
{
  myRoot = NULL;
  myTolerance = tol;
  maxChildren = 8;
  minChildren = 2;
}
template <class Z> MY_INLINE RTree<Z>::RTree (double tol, int max_c, int min_c)
{
  myRoot = NULL;
  myTolerance = tol;
  maxChildren = max_c;
  minChildren = min_c;
}
template <class Z> MY_INLINE RTree<Z>::~RTree()
{
  if ( myRoot != NULL )
  {
      //Go through and get all the children in a list.
    DLIList <RTreeNode<Z>*> to_delete;
    to_list(to_delete, myRoot);
    int ii;
    for(ii = to_delete.size(); ii > 0; ii-- )
      delete to_delete.pop();
    delete myRoot;
  }
}
template <class Z> MY_INLINE void RTree<Z>::to_list(DLIList <RTreeNode<Z>*> &member_list,
                                          RTreeNode<Z> *top)
{
    //Get the children of the top into the list.
  int ii;
  RTreeNode <Z> *curr_node;
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

template <class Z> MY_INLINE CubitStatus RTree<Z>::add(Z data)
{
  CubitStatus stat;
  RTreeNode<Z> *new_root;
  RTreeNode<Z> *new_node = new RTreeNode<Z> (data, myTolerance, maxChildren,
                                             minChildren);
  if ( myRoot == NULL )
  {
    CubitBox b_box = data->bounding_box();
    myRoot = new RTreeNode <Z> (b_box, maxChildren, minChildren);
	myRoot->set_leaf_level(LEAF_RNODE);
    stat = myRoot->insert(new_node, new_root);
      //this shouldn't change the root, or fail!
    if ( stat != CUBIT_SUCCESS || new_root != NULL )
    {
      PRINT_ERROR("Insertion into RTree failed.\n");
      return CUBIT_FAILURE;
    }
    return CUBIT_SUCCESS;
  }
  
  stat = myRoot->insert(new_node, new_root);
  if ( stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Insertion into RTree failed.\n");
    return CUBIT_FAILURE;
  }
  if ( new_root != NULL )
  {
      //this is fine, it just means we are adding more
      //so the root had to be split...
    myRoot = new_root;
  }
  return CUBIT_SUCCESS;
}
template <class Z> MY_INLINE CubitStatus RTree<Z>::find(const CubitBox &range_box,
                                              DLIList <Z> &range_members )
{
    //Find all of the members of the RTree that intersect this range_box.
  if ( myRoot == NULL )
  {
      // Nothing has been added to this Tree yet, so we are not going to find this
      // object in it.
      return CUBIT_SUCCESS;
  }
  CubitStatus stat = recursive_find(myRoot, range_box, range_members);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  else
    return CUBIT_SUCCESS;
}
template <class Z> MY_INLINE CubitStatus RTree<Z>::recursive_find(RTreeNode<Z> *rect_tree,
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
  RTreeNode<Z> *curr_node;
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
template <class Z> MY_INLINE CubitBoolean RTree<Z>::remove( Z data )
{
  RTreeNode<Z> *new_root = NULL;
  CubitBoolean delete_root = CUBIT_FALSE;
  CubitBoolean return_val = myRoot->remove( data, new_root, delete_root );
  if ( new_root != NULL )
  {
      //Only if we are condensing the tree do we want to delete the root node.
      //There are other reasons the root has changed (rebalance...), in which
      //cases the root is now a child of the new root...
    if ( delete_root )
      delete myRoot;
    myRoot = new_root;
  }
  return return_val;
}
//--------------------------------------------------------------------------
//Algorithm: min_dist_sq
//Description:  Finds the minimum distance squared between the given
//              point and the box. If the point is on or in the box, the
//              min distance is zero.
//--------------------------------------------------------------------------
template <class Z> MY_INLINE
double RTree<Z>::min_dist_sq(CubitVector &q,
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
template <class Z> MY_INLINE
bool RTree<Z>::less_than_func(RTreeNode<Z> *&node_a,
                              RTreeNode<Z> *&node_b)
{
  if ( node_a->get_dist() < node_b->get_dist() )
    return true;
  else
    return false;
}

template <class Z> MY_INLINE
CubitStatus RTree<Z>::k_nearest_neighbor(CubitVector &q,
                                         int k,
                                         double &closest_dist,
                                         DLIList<Z> &nearest_neighbors,
                                         typename RTree<Z>::DistSqFunc dist_sq_point_data)
{
    //first create the priority queue.
  PriorityQueue< RTreeNode<Z>*> near_queue(RTree<Z>::less_than_func);
  
  myRoot->set_dist(0.0);
  near_queue.push(myRoot);
  RTreeNode<Z> *element, *child_element;
  int num_found = 0;
  int ii;
  double data_dist, box_dist;
  Z data;
  
  while( !near_queue.empty() )
  {
    element = near_queue.top();
    near_queue.pop();
    if ( element->is_data() )
    {
      data = element->get_data();
        //calculate the exact distance.
      data_dist = dist_sq_point_data(q, data);
        //compare this distance with the next item's distance.
      if ( element->dist_is_box() && !near_queue.empty() &&
           data_dist > near_queue.top()->get_dist())
      {
          //If its bigger, add it back into the list
          //with the updated distance.
        element->set_dist(data_dist);
        near_queue.push(element);
        element->set_dist_is_box(0);
      }
      else
      {
        nearest_neighbors.append(element->get_data());
        if ( num_found == 0 )
          closest_dist = element->get_dist();
        num_found++;
        if ( num_found == k )
          return CUBIT_SUCCESS;
      }
    }
    else 
    {
      for ( ii = 0; ii < element->num_children(); ii++ )
      {
        child_element = element->get_child(ii);
        CubitBox bounding_box = child_element->bounding_box();
        box_dist = min_dist_sq(q, bounding_box);
        child_element->set_dist(box_dist);
        near_queue.push(child_element);
      }
    }
  }
  return CUBIT_FAILURE;
}

    
    
